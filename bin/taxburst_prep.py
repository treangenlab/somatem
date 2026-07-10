#!/usr/bin/env python3

# Author: Austin Marshall
# Date: 12/Sep/25


import argparse, csv, json, math
from collections import defaultdict

RANKS_IN = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
LEMUR_COLS = {
    "species": "species",
    "genus": "genus",
    "family": "family",
    "order": "order",
    "class": "class",
    "phylum": "phylum",
    # "clade" is present in LEMUR but not needed for Krona/Taxburst
    "superkingdom": "superkingdom",
    "fraction": r"^(F|abundance)$"  # Match either 'F' or 'abundance' for fraction
}

def sanitize_taxonomy_for_taxburst(out):
    """
    Ensure Taxburst/Krona never sees empty taxonomy names.
    Fills blanks with deterministic placeholders using parent context.
    """
    # Strip leading/trailing whitespace from every taxonomy rank value so that
    # blank-padded fields (common in CSV exports) are treated as truly empty
    for k in RANKS_IN:
        v = out.get(k, "")
        out[k] = v.strip() if isinstance(v, str) else ""

    # Top level must never be empty
    if not out["superkingdom"]:
        out["superkingdom"] = "unclassified"

    # Fill any remaining missing ranks with non-empty placeholders
    # Include parent context to reduce risk of duplicate-name assertions in Taxburst
    for i, rank in enumerate(RANKS_IN[1:], start=1):
        if not out[rank]:
            parent_rank = RANKS_IN[i - 1]
            parent_name = out[parent_rank] or "root"
            # keep placeholder compact and filesystem/html safe-ish
            parent_safe = str(parent_name).replace(" ", "_")
            out[rank] = f"unclassified_{rank}_in_{parent_safe}"

    return out

def detect_delimiter(path):
    # Default to tab; fall back to comma if needed
    with open(path, "r", newline="") as fh:
        head = fh.read(4096)
    return "\t" if "\t" in head and head.count("\t") >= head.count(",") else ","

def read_lemur_rows(path, from_tool):
    delim = detect_delimiter(path)
    with open(path, "r", newline="") as fh:
        r = csv.DictReader(fh, delimiter=delim)
        # normalize header keys (strip spaces)
        fieldmap = {k.strip(): k for k in r.fieldnames}
        def get(row, key):
            src = LEMUR_COLS.get(key, key)
            # If the source is a regex pattern (starts and ends with /)
            if key == "fraction" and src.startswith("^") and src.endswith("$"):
                import re
                pattern = re.compile(src)
                for field in fieldmap.values():
                    if pattern.match(field):
                        return (row.get(field, "") or "").strip()
                return ""
            # Normal dictionary lookup
            src = fieldmap.get(src, src)
            return (row.get(src, "") or "").strip()

        for row in r:
            # Get the tax_id value
            tax_id = get(row, 'tax_id')

            # Skip rows with 'unmapped' or 'unclassified' in tax_id (for EMU)
            if from_tool == "emu":
                if tax_id and ('unmapped' in tax_id.lower() or 'unclassified' in tax_id.lower()):
                    continue
            
            out = {}
            for k in RANKS_IN:
                out[k] = get(row, k)

            # replace class of `Actinobacteria` to `Actinomycetes` in the row (for EMU) : prevent duplicate with phylyum
            if from_tool == "emu":
                if out["class"] == "Actinobacteria":
                    out["class"] = "Actinomycetes"

            out = sanitize_taxonomy_for_taxburst(out)

            # parse fraction (abundance) from LEMUR 'F'
            f_str = get(row, "fraction")
            try:
                frac = float(f_str)
            except Exception:
                frac = 0.0
            out["fraction"] = max(0.0, frac)
            yield out

def read_sylph_rows(path):
    """Convert leaf clades from sylph-tax's .sylphmpa format."""
    with open(path, "r", newline="") as fh:
        records = list(csv.DictReader(
            (line for line in fh if not line.startswith("#")),
            delimiter="\t"
        ))

    clades = [(row.get("clade_name") or "").strip() for row in records]
    leaf_clades = {
        clade for clade in clades if clade and not any(
            other.startswith(clade + "|") for other in clades
        )
    }
    prefix_to_rank = {
        "d": "superkingdom",
        "r": "superkingdom",
        "k": "superkingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }
    aggregated = defaultdict(float)

    for row in records:
        clade = (row.get("clade_name") or "").strip()
        if clade not in leaf_clades:
            continue

        taxonomy = {rank: "" for rank in RANKS_IN}
        for component in clade.split("|"):
            if "__" not in component:
                continue
            prefix, name = component.split("__", 1)
            rank = prefix_to_rank.get(prefix)
            if rank and name:
                taxonomy[rank] = name

        taxonomy = sanitize_taxonomy_for_taxburst(taxonomy)
        try:
            fraction = float(row.get("relative_abundance") or 0.0) / 100.0
        except (TypeError, ValueError):
            fraction = 0.0
        aggregated[tuple(taxonomy[rank] for rank in RANKS_IN)] += max(0.0, fraction)

    for lineage, fraction in aggregated.items():
        out = dict(zip(RANKS_IN, lineage))
        out["fraction"] = fraction
        yield out

def write_krona(rows, out_path):
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        # Krona-style header that taxburst accepts
        w.writerow(["fraction"] + RANKS_IN)
        for r in rows:
            w.writerow([f"{r['fraction']:.12g}"] + [r[k] for k in RANKS_IN])

def rows_to_json_nodes(rows, multiplier=1000.0):
    # Build a tree keyed by (rank, name)
    def make_node(name, rank):
        return {"name": name, "rank": rank, "count": 0.0, "children": {}}

    root_by_name = {}  # for top-level ranks (superkingdom)

    for r in rows:
        count = r["fraction"] * multiplier
        if count <= 0:
            continue
        parent_map = root_by_name
        parent_node = None

        # Walk the taxonomy from superkingdom -> ... -> species
        for i, rank in enumerate(RANKS_IN):
            name = r[rank]
            if not name:
                # skip missing rank; continue to next rank
                continue
            # choose where to store this node
            store = parent_map
            key = (rank, name)
            if key not in store:
                node = make_node(name, rank.capitalize())
                store[key] = node
            else:
                node = store[key]
            # accumulate count at every node along the path
            node["count"] += count
            # descend
            parent_node = node
            parent_map = node["children"]

    # convert dict-children to lists recursively
    def finalize(node_dict):
        out = []
        for (_rk, _nm), nd in node_dict.items():
            nd_final = {
                "name": nd["name"],
                "rank": nd["rank"],
                # round to integer counts like other parsers (1000 × fraction)
                "count": int(round(nd["count"]))
            }
            kids = finalize(nd["children"])
            if kids:
                nd_final["children"] = kids
            out.append(nd_final)
        return out

    return finalize(root_by_name)

def main():
    ap = argparse.ArgumentParser(
        description="Convert LEMUR taxonomy TSV to Taxburst input (Krona TSV or JSON)."
    )
    ap.add_argument("input", help="Taxonomic profile generated by the selected tool")
    ap.add_argument("--from_tool", choices=["lemur", "emu", "sylph"], default="lemur", help="The tool generating the taxonomy input")
    ap.add_argument("-o", "--output", required=True, help="Output file path")
    ap.add_argument("-F", "--format", choices=["krona", "json"], default="krona",
                    help="Taxburst input format to write")
    ap.add_argument("--multiplier", type=float, default=1000.0,
                    help="Counts = fraction × multiplier for JSON output (default 1000)")
    args = ap.parse_args()

    rows = list(read_sylph_rows(args.input) if args.from_tool == "sylph" else read_lemur_rows(args.input, args.from_tool))

    if args.format == "krona":
        write_krona(rows, args.output)
    else:
        nodes = rows_to_json_nodes(rows, multiplier=args.multiplier)
        with open(args.output, "w") as fh:
            json.dump(nodes, fh, indent=2)

if __name__ == "__main__":
    main()
