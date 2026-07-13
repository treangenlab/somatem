#!/usr/bin/env python3
"""
Generates html output which shows threat potential for each sequence.
Also shows all information about each sequence.
"""
import argparse
import os
import sys
import time
import json
import csv
import re
from html import escape
from collections import defaultdict
from collections import Counter
from distutils.dir_util import copy_tree
try:
    from jinja2 import Environment, FileSystemLoader
except ImportError:
    Environment = None
    FileSystemLoader = None
# Updated the cvs field size to resolve "_csv.Error: field larger than field limit"
goTerms = {}
funsocs_start_idx = 5
funsocs_end_idx = 37
funsocs_col_names = {}

num_cols = 7
mode = 'Default'
# ===========================================================
# Classes

class Sequence:
    """
    Information that can be extracted for each single sequence.
    """

    def __init__(self, query, length, tax_id, uniprot, gene_name, sequence, go_childs, go_terms, funsocs_dict,
            hsp=None, aligned=False, bsat_hit=False, vfdb_hit="NA", flag="0"):
        self.length = length
        self.query = query
        self.tax_id = tax_id
        self.uniprot = uniprot
        self.gene_name = gene_name
        self.sequence = sequence
        self.go_childs = go_childs
        self.go_terms = go_terms
        self.hsp = hsp
        self.aligned = aligned
        self.bsat_hit = bsat_hit
        self.vfdb_hit = vfdb_hit
        self.flag = flag
        self.funsocs = {}
        for val in funsocs_dict:
            self.funsocs[val] = funsocs_dict[val]

class HSPhit:
    """
    Information provided by blast output.
    """

    def __init__(self):
        self.qSeq = ""
        self.hSeq = ""
        self.midLine = ""
        self.bit_score = 0
        self.evalue = 100000
        self.hitFrom = 0
        self.hitTo = 0
        self.queryFrom = 0
        self.queryTo = 0
        self.hit = ""


# ===========================================================
# OS operations

def create_folder(path):
    """
    Creates a folder at a given path.
    :param path: the path where the folder should be created. Includes the path and the folder name
    """
    if not os.path.exists(path):
        os.makedirs(path)


def create_output_file(t_file, o_file, output_prefix, funsocs_count, version, funsocs_names):
    """
    Creates the output file by using a html template and filling it with the necessary values.
    :param t_file: path to the html template file
    :param o_file: path where the output file should be created
    :param output_prefix: the prefix for the output
    :param funsocs_count: a dictionary counting the occurences of each funsocs
    :param version: string representation of the version that is used
    :param mtimes: list with last-modified times of input files
    :param funsocs_names: the names of the funsocs column
    """
    global num_cols

    # create jinja2 environment to write output file with the html template and calculated vars.
    # env = Environment(loader=FileSystemLoader(os.getcwd()))
    folder, template = os.path.split(t_file)
    env = Environment(loader=FileSystemLoader(folder))

    # assign variables to defined variables in the html template file
    # result = env.get_template(os.path.relpath(t_file)).render(
    result = env.get_template(template).render(
        output_prefix=f'"{output_prefix}"',
        funsocs_count=funsocs_count,
        version=version,
        funsocs_names=funsocs_names,
        column_count=funsocs_end_idx-funsocs_start_idx + num_cols
    )

    # create the output file
    with open(o_file, 'w') as file:
        file.write(result)


def create_executive_report(o_file, sequences, funsocs_count, funsocs_names, version, mode, rflag, fasta_file, mtimes):
    """
    Create a plain-language landing page for people who need the take-home result
    before digging into the technical table.
    """
    total = len(sequences)
    concern_sequences = []
    organism_counts = Counter()
    gene_counts = Counter()
    function_counts = Counter()
    function_by_taxon = defaultdict(Counter)
    safety_reason_counts = Counter()

    def truthy(value):
        return str(value).strip().lower() in ["yes", "y", "1", "true", "flag", "flagged"]

    def clean_label(value, empty="Unassigned"):
        value = str(value or "").strip()
        return value if value and value != "-" else empty

    def display_function(name):
        name = str(name or "").replace("_", " ").strip()
        return re.sub(r"(\w)([A-Z])", r"\1 \2", name).title()

    def sequence_functions(seq):
        found = []
        for idx in sorted(funsocs_col_names):
            if str(seq.funsocs.get(idx, "0")).strip() == "1":
                found.append(display_function(funsocs_col_names[idx]))
        return found

    function_notes = {}
    for item in funsocs_names:
        function_notes[display_function(item.get("name"))] = item.get("title")

    for seq in sequences:
        taxon = clean_label(seq.tax_id)
        functions = sequence_functions(seq)
        reasons = []

        if truthy(seq.bsat_hit):
            reasons.append("BSAT/select-agent database hit")
        if truthy(seq.vfdb_hit):
            reasons.append("Virulence factor database hit")
        if truthy(seq.flag):
            reasons.append("SeqScreen flagged this sequence")
        if functions:
            reasons.append("Functional risk signal")

        if reasons:
            concern_sequences.append((seq, taxon, functions, reasons))
            safety_reason_counts.update(reasons)

        organism_counts[taxon] += 1
        if seq.gene_name and seq.gene_name != "-":
            gene_counts[seq.gene_name] += 1
        for function in functions:
            function_counts[function] += 1
            function_by_taxon[taxon][function] += 1

    concern_count = len(concern_sequences)
    concern_percent = round((concern_count / total) * 100, 1) if total else 0
    blacklist_count = safety_reason_counts["BSAT/select-agent database hit"] + safety_reason_counts["Virulence factor database hit"]
    top_organisms = organism_counts.most_common(8)
    top_genes = gene_counts.most_common(8)

    if concern_count == 0:
        status_label = "No immediate safety concern signals detected"
        status_class = "ok"
        status_text = "SeqScreen did not find blacklist, virulence-factor, or functional risk signals in this sample."
    elif blacklist_count > 0:
        status_label = "Review recommended before release or sharing"
        status_class = "alert"
        status_text = "At least one sequence matched a safety-focused database. Treat this as a screening signal and have a qualified reviewer inspect the technical results."
    else:
        status_label = "Functional signals detected"
        status_class = "watch"
        status_text = "Some sequences had functions associated with biological activity. These are not diagnoses, but they are worth reviewing in context."

    def list_items(items, empty_text):
        if not items:
            return f"<li>{escape(empty_text)}</li>"
        return "".join(
            f"<li><span>{escape(str(name))}</span><strong>{count}</strong></li>"
            for name, count in items
        )

    function_rows = "".join(
        f"<tr><td>{escape(function)}</td><td>{count}</td><td>{escape(str(function_notes.get(function) or ''))}</td></tr>"
        for function, count in function_counts.most_common()
    ) or "<tr><td colspan=\"3\">No functional signals were detected.</td></tr>"

    function_taxon_rows = []
    for taxon, counts in sorted(function_by_taxon.items(), key=lambda item: sum(item[1].values()), reverse=True):
        for function, count in counts.most_common():
            function_taxon_rows.append(
                f"<tr><td>{escape(taxon)}</td><td>{escape(function)}</td><td>{count}</td></tr>"
            )
    function_taxon_table = "".join(function_taxon_rows[:80]) or "<tr><td colspan=\"3\">No function-by-taxon signals were detected.</td></tr>"

    review_rows = []
    for seq, taxon, functions, reasons in concern_sequences[:40]:
        review_rows.append(
            "<tr>"
            f"<td>{escape(str(seq.query))}</td>"
            f"<td>{escape(taxon)}</td>"
            f"<td>{escape(clean_label(seq.gene_name, 'No gene name'))}</td>"
            f"<td>{escape(', '.join(functions) if functions else 'No functional category')}</td>"
            f"<td>{escape('; '.join(reasons))}</td>"
            "</tr>"
        )
    review_table = "".join(review_rows) or "<tr><td colspan=\"5\">No sequences require review based on these screening signals.</td></tr>"

    input_rows = "".join(
        f"<tr><td>{escape(str(name))}</td><td>{escape(str(modified))}</td></tr>"
        for name, modified in mtimes.items()
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>SeqScreen Summary</title>
  <style>
    :root {{
      --ink: #18212f;
      --muted: #5c6675;
      --line: #d9dee7;
      --bg: #f6f8fb;
      --panel: #ffffff;
      --ok: #0f766e;
      --watch: #9a6700;
      --alert: #b42318;
      --blue: #235789;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--ink);
      background: var(--bg);
      line-height: 1.5;
    }}
    header {{
      background: #ffffff;
      border-bottom: 1px solid var(--line);
      padding: 28px 32px 20px;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 24px 20px 48px;
    }}
    h1, h2, h3 {{ margin: 0; }}
    h1 {{ font-size: 2rem; }}
    h2 {{ font-size: 1.15rem; margin-bottom: 12px; }}
    p {{ margin: 8px 0 0; }}
    .subtle {{ color: var(--muted); }}
    .status {{
      margin-top: 18px;
      padding: 18px;
      border: 1px solid var(--line);
      background: var(--panel);
      border-left: 8px solid var(--blue);
    }}
    .status.ok {{ border-left-color: var(--ok); }}
    .status.watch {{ border-left-color: var(--watch); }}
    .status.alert {{ border-left-color: var(--alert); }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 14px;
      margin: 20px 0;
    }}
    .metric, section {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 18px;
    }}
    .metric strong {{
      display: block;
      font-size: 2rem;
      line-height: 1;
      margin-bottom: 8px;
    }}
    .two-col {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 16px;
      margin-bottom: 16px;
    }}
    .wide {{ margin-bottom: 16px; }}
    ul.clean {{
      list-style: none;
      padding: 0;
      margin: 0;
    }}
    ul.clean li {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      padding: 9px 0;
      border-top: 1px solid #edf0f5;
    }}
    ul.clean li:first-child {{ border-top: 0; }}
    a.button {{
      display: inline-block;
      margin: 12px 10px 0 0;
      padding: 10px 14px;
      color: #fff;
      background: var(--blue);
      border-radius: 6px;
      text-decoration: none;
      font-weight: 700;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.92rem;
      background: #fff;
    }}
    th, td {{
      border-top: 1px solid #edf0f5;
      padding: 8px;
      vertical-align: top;
      text-align: left;
    }}
    th {{
      background: #eef3f8;
      color: #2b3443;
      font-weight: 700;
    }}
    .note {{
      color: var(--muted);
      font-size: 0.94rem;
    }}
    .badge {{
      display: inline-block;
      padding: 4px 8px;
      border-radius: 999px;
      background: #eef3f8;
      color: #2b3443;
      font-weight: 700;
      font-size: 0.82rem;
    }}
    @media (max-width: 860px) {{
      .grid, .two-col {{ grid-template-columns: 1fr; }}
      header {{ padding: 22px 20px; }}
      table {{ display: block; overflow-x: auto; white-space: nowrap; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>SeqScreen Summary</h1>
    <p class="subtle">Plain-language overview of pathogen and functional-risk screening results. This file is self-contained and can be shared by itself.</p>
    <div class="status {status_class}">
      <h2>{escape(status_label)}</h2>
      <p>{escape(status_text)}</p>
    </div>
  </header>
  <main>
    <div class="grid">
      <div class="metric"><strong>{total}</strong><span>Total sequences screened</span></div>
      <div class="metric"><strong>{concern_count}</strong><span>Sequences with concern signals</span></div>
      <div class="metric"><strong>{concern_percent}%</strong><span>Share of screened sequences flagged</span></div>
      <div class="metric"><strong>{blacklist_count}</strong><span>Safety database matches</span></div>
    </div>

    <section>
      <h2>How To Read This</h2>
      <p>SeqScreen looks for signals associated with known pathogens, virulence, antimicrobial resistance, and related biological functions. A flagged result is not a final diagnosis. It means the sequence should be reviewed by a domain expert before operational decisions are made.</p>
      <p class="note">Plain-language rule of thumb: safety database matches are the highest-priority review items; functional signals explain what kind of biological activity was detected; taxa show which organism assignments those signals were associated with.</p>
    </section>

    <div class="two-col" style="margin-top:16px;">
      <section>
        <h2>Most Common Organisms</h2>
        <ul class="clean">{list_items(top_organisms, "No organism assignments were reported.")}</ul>
      </section>
      <section>
        <h2>Most Common Gene Names</h2>
        <ul class="clean">{list_items(top_genes, "No gene names were reported.")}</ul>
      </section>
    </div>

    <section class="wide">
      <h2>Functions Found</h2>
      <p class="note">These are the functional categories detected by SeqScreen, sorted by how often they appeared.</p>
      <table><thead><tr><th>Function</th><th>Sequences</th><th>Plain-language note</th></tr></thead><tbody>{function_rows}</tbody></table>
    </section>

    <section class="wide">
      <h2>Functions By Taxon</h2>
      <p class="note">This table answers: what functions were found, and which organism assignments were they associated with?</p>
      <table><thead><tr><th>Taxon</th><th>Function</th><th>Sequences</th></tr></thead><tbody>{function_taxon_table}</tbody></table>
    </section>

    <section class="wide">
      <h2>Sequences To Review</h2>
      <p class="note">This is a short triage list for reviewers. It is capped at the first 40 review items to keep the shared report readable.</p>
      <table><thead><tr><th>Sequence</th><th>Taxon</th><th>Gene</th><th>Functions</th><th>Why It Is Listed</th></tr></thead><tbody>{review_table}</tbody></table>
    </section>

    <section class="wide">
      <h2>Run Details</h2>
      <p><span class="badge">SeqScreen {escape(version)}</span></p>
      <p><strong>Mode:</strong> {escape(mode)}</p>
      <p><strong>Runtime flags:</strong> {escape(rflag)}</p>
      <p><strong>Input FASTA:</strong> {escape(os.path.basename(fasta_file))}</p>
    </section>

    <section>
      <h2>Input Files</h2>
      <table>{input_rows}</table>
    </section>
  </main>
</body>
</html>
"""
    with open(o_file, 'w') as file:
        file.write(html)


def create_go_file(template, outfile, sequences, go_network):
    """
    Create output for go term visualization from template and save to outfile.
    :param template: path to html template file
    :param outfile: output file path
    :param sequences: formatted sequences
    :param go_network: list of {from: go_id, to: go_id} directed graph edges
    """
    # create jinja2 environment
    folder, template = os.path.split(template)
    env = Environment(loader=FileSystemLoader(folder))

    # pass variables to template
    result = env.get_template(template).render(
        sequences=sequences,
        go_network=go_network
    )

    # write output file
    with open(outfile, 'w') as file:
        file.write(result)


def create_fallback_output_file(o_file, funsocs_count, version, funsocs_names):
    """
    Create a small technical report if the legacy SeqScreen template is missing.
    """
    rows = "".join(
        f"<tr><td>{escape(str(item['name']).replace('_', ' '))}</td><td>{escape(str(item['count']))}</td></tr>"
        for item in funsocs_count
    ) or "<tr><td colspan=\"2\">No functional concern categories were detected.</td></tr>"
    descriptions = "".join(
        f"<li><strong>{escape(str(item['name']))}</strong>: {escape(str(item['title']))}</li>"
        for item in funsocs_names
        if item.get('title')
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>SeqScreen Technical Report</title>
  <style>
    body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; color: #18212f; background: #f6f8fb; line-height: 1.5; }}
    main {{ max-width: 980px; margin: 0 auto; padding: 28px 20px 44px; }}
    h1 {{ margin: 0 0 6px; font-size: 2rem; }}
    h2 {{ margin-top: 24px; font-size: 1.2rem; }}
    a {{ color: #235789; font-weight: 700; }}
    table {{ width: 100%; border-collapse: collapse; background: #fff; border: 1px solid #d9dee7; }}
    td, th {{ padding: 10px; border-top: 1px solid #edf0f5; text-align: left; }}
    th {{ background: #eef3f8; }}
  </style>
</head>
<body>
  <main>
    <h1>SeqScreen Technical Report</h1>
    <p>SeqScreen version: {escape(str(version))}</p>
    <p>The legacy technical template was not found, so this compact report was generated instead. The full parsed data are still available in the <code>data/</code> folder.</p>
    <p><a href="index.html">Back to plain-language summary</a> | <a href="go.html">Open GO term view</a></p>
    <h2>Functional Concern Categories</h2>
    <table><thead><tr><th>Category</th><th>Sequences</th></tr></thead><tbody>{rows}</tbody></table>
    <h2>Category Notes</h2>
    <ul>{descriptions or "<li>No category descriptions were available.</li>"}</ul>
  </main>
</body>
</html>
"""
    with open(o_file, 'w') as file:
        file.write(html)


def create_fallback_go_file(outfile, go_network):
    """
    Create a minimal GO report if the legacy GO visualization template is missing.
    """
    edge_count = len(go_network) if go_network else 0
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>SeqScreen GO Terms</title>
  <style>
    body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; color: #18212f; background: #f6f8fb; line-height: 1.5; }}
    main {{ max-width: 820px; margin: 0 auto; padding: 28px 20px 44px; }}
    a {{ color: #235789; font-weight: 700; }}
  </style>
</head>
<body>
  <main>
    <h1>GO Term View</h1>
    <p>The legacy GO visualization template was not found, so the interactive graph could not be rendered.</p>
    <p>{edge_count} GO network relationships were loaded for this run.</p>
    <p><a href="index.html">Back to plain-language summary</a> | <a href="out.html">Open technical report</a></p>
  </main>
</body>
</html>
"""
    with open(outfile, 'w') as file:
        file.write(html)


def write_file(sequences, archive_name, mode, rflag, fasta_file, length,mtimes):
    
    query_id = 0
    data = []
    count = 1  
    summary = {}  
    summary['mode'] = mode
    summary['rflag'] = rflag
    summary['fasta_file']= fasta_file
    summary['length'] = length
    summary['mtimes'] = mtimes
    data.append(summary)
    # Specify the folder for the split output json files
    out_json_folder = os.path.abspath(archive_name)+'/data'
    create_folder(out_json_folder)

    for seq in sequences:

        item = {}
        item['length'] = seq.length
        item['query'] = seq.query
        item['taxID'] = seq.tax_id
        item['uniprot'] = seq.uniprot
        item['geneName'] = seq.gene_name
        item['sequence'] = seq.sequence
        item['qID'] = query_id
        item['go_terms'] = seq.go_terms
        item['go_childs'] = seq.go_childs
        item['hasAlignment'] = str(seq.aligned).lower()
        item['bsat'] = seq.bsat_hit
        item['vfdb'] = seq.vfdb_hit
        item['flag'] = seq.flag
        range_val = funsocs_end_idx - funsocs_start_idx +1 
        for val in range(range_val):
            item[funsocs_col_names[val]] = seq.funsocs[val] 

        if seq.aligned:
            item['alignment'] = vars(seq.hsp)

        query_id += 1
        data.append(item)
        
        if(query_id % 25000 == 0):
            data[0]['length']= 25000
            file_name = "output_part"+str(count)+".json"
            seq_file = os.path.join(out_json_folder, file_name)
            with open(seq_file,'w')as seqfile:
                json.dump(data,seqfile)
            data = []
            data.append(summary)
            count += 1

    if(query_id % 25000 != 0):
        file_name = "output_part"+str(count)+".json"
        data[0]['length']= (query_id % 25000)
        seq_file = os.path.join(out_json_folder, file_name)
        with open(seq_file,'w')as seqfile:
            json.dump(data,seqfile)
        data = []
        data.append(summary)
def format_go_information(sequences):
    """
    Formats the sequence information leaving only relevant information.
    :param sequences: list/array of sequences
    :return: the formatted list of sequences.
    """
    data = []
    query_id = 0

    for seq in sequences:
        item = {}
        item['query'] = seq.query
        item['qID'] = query_id
        item['go_terms'] = seq.go_terms
        item['go_childs'] = seq.go_childs
        query_id += 1
        data.append(item)

    return data

def format_sequence_information(sequences):
    """
    Formats the sequence information leaving only relevant information.
    :param sequences: list/array of sequences
    :return: the formatted list of sequences.
    """
    data = []
    query_id = 0

    for seq in sequences:

        item = {}
        item['query'] = seq.query
        item['length'] = seq.length
        item['taxID'] = seq.tax_id
        item['uniprot'] = seq.uniprot
        item['geneName'] = seq.gene_name
        item['sequence'] = seq.sequence
        item['qID'] = query_id
        item['go_terms'] = seq.go_terms
        item['go_childs'] = seq.go_childs
        item['hasAlignment'] = str(seq.aligned).lower()
        item['bsat'] = seq.bsat_hit
        range_val = funsocs_end_idx - funsocs_start_idx +1 
        for val in range(range_val):
            item[funsocs_col_names[val]] = seq.funsocs[val] 

        if seq.aligned:
            item['alignment'] = vars(seq.hsp)
        query_id += 1
        data.append(item)

    return data


# ===========================================================
# Count funsocss
def build_funsocs_count(tally):
    ret = []
    for funsocs in tally:
        ret.append({
            'name': funsocs,
            'count': tally[funsocs],
            'style': 'light'
        })
    return ret

def build_funsocs_names(tally,title):
    ret = []
    for funsocs in tally:
        name = funsocs['name']
        name = name.replace("_"," ")
        name = re.sub(r"(\w)([A-Z])", r"\1 \2", name)
        ret.append({
              'name':name,
               'title':title[funsocs['name']]
        })
    return ret

def count_funsocs(sequences):
    tally = defaultdict(int)
    val = funsocs_end_idx - funsocs_start_idx+1
    for seq in sequences:
        for idx in range(val):
            tally[funsocs_col_names[idx]] += seq.funsocs[idx] == '1'
    return tally

def set_title(funsocs_file):
    if not funsocs_file:
        return defaultdict(int)
    ffile = csv.DictReader(open(funsocs_file, 'r'), delimiter='\t')
    title = defaultdict(int)
    for line in ffile:
        idx = line["Name"]
        title[idx]= line["Description"]
    return title

# ===========================================================
# Load stuff

def load_all(fasta):
    """
    Loads all sequences with their names and length into a dictionary
    :param fasta: path to the fasta file
    :return: dictionary with sequence names and lengths
    """
    seq_names_lengths = {}
    last = ""
    with open(fasta, 'r') as file:
        for line in file:
            # rstrip removes whitespace before and after the given string
            line = line.strip()
            if line.startswith('>'):
                last = line[1:].split(' ')[0]
                seq_names_lengths[last] = ""
            elif last in seq_names_lengths:
                seq_names_lengths[last] += line
    return seq_names_lengths


def load_sequences(report_file, seq_names_lengths, alignments, go_names):
    """
    Combines all extracted information to get a list with all sequences and their information
    :param report_file: output file from the SeqScreen pipeline
    :param seq_names_lengths: dictionary with name and length of each sequence
    :param alignments: alignments extracted from the blast output file
    :param go_names: dict with 'go-number': go-name - elements
    :return: list with all sequences and their specific information
    """
    ret = []
    visited = set()
    funsocs = {}
    inputs = open(report_file, 'r').readlines()
    field_map = define_column_positions(inputs[0])
    
    query_idx = field_map["query"]
    size_idx = field_map["size"]
    organism_idx = field_map["organism"]
    gene_idx = field_map["gene_name"]
    uniprot_idx = field_map["uniprot"]
    go_idx = field_map["go"]
    bsat_idx = field_map.get("bsat_hit")
    vfdb_idx = field_map.get("vfdb_hit")
    flag_idx = field_map.get("flag")

    for line in inputs[1:]:
        fields = line.rstrip("\n").split("\t")
        # GO Terms: create list of {'id', 'name'} elements from fields["go_idx"] and go_names dict
        go_childs = list(map(lambda term: {
            'id': term.strip(),
            'name': go_names.get(term.strip(), 
                'No molecular function or biological process Go terms found')
        },fields[go_idx].split(";")))

        goTerms = getAllGoTerms(fields[go_idx].split(";"))
        fields[go_idx] = list(map(lambda term: {
            'id': term.strip(),
            'name': go_names.get(term.strip(), 
                'No molecular function or biological process Go terms found')
        }, goTerms))

        bsat = fields[bsat_idx] if bsat_idx is not None else 'NA'
        vfdb = fields[vfdb_idx] if vfdb_idx is not None else 'NA'
        flag = fields[flag_idx] if flag_idx is not None else '0'
       
        # Updating value of funsocs column in the funsocs array
        idx = 0
        for index in range(funsocs_start_idx, funsocs_end_idx+1):
            funsocs[idx] = fields[index] if index < len(fields) else "0"
            idx = idx+1
       
        if fields[query_idx] in alignments:
            hsp = alignments[fields[query_idx]]
            ret.append(Sequence(fields[query_idx], fields[size_idx],
                                fields[organism_idx], fields[uniprot_idx], fields[gene_idx], 
                                seq_names_lengths[fields[query_idx]], go_childs, fields[go_idx], funsocs, hsp, True, bsat, vfdb, flag))
        else:
            ret.append(Sequence(fields[query_idx], fields[size_idx],
                                fields[organism_idx], fields[uniprot_idx], fields[gene_idx],
                                seq_names_lengths[fields[query_idx]], go_childs, fields[go_idx], funsocs, bsat_hit=bsat, vfdb_hit=vfdb, flag=flag))

        visited.add(fields[query_idx])

    for gene in list(set(seq_names_lengths) - visited):
        ret.append(Sequence(gene, len(seq_names_lengths[gene]), "-", "-", "-",
                            seq_names_lengths[gene], None, None, defaultdict(int)))

    return ret


def load_go_names(names_file):
    """
    Load all go terms with their names into a dict
    :param names_file: path to the go_names file
    :return: dict with idx: name
    """
    go_names = {}
    if not names_file:
        return go_names
    with open(names_file, 'r') as file:
        for line in file:
            line = line.strip()
            idx, name = line.split("\t")

            if id not in go_names:
                go_names[idx] = name

    return go_names


def load_go_network(network_file):
    """
    Load go network into a list
    :param network_file: path to the go_network file
    :return: list with {from: id, to: id} elements
    """
    global goTerms
    go_network = []
    if not network_file:
        return go_network
    with open(network_file, 'r') as file:
        for line in file:
            line = line.strip()
            child, parent = line.split("\t")

            go_network.append({'from': child, 'to': parent, 'arrows': 'from'})
            if(parent in goTerms):
                goTerms[parent].append(child)
            else:
                goTerms[parent] = [child]

    return go_network


# ===========================================================
# Alignments

def extract_alignments(xml_file):
    """
    Extracts the information provided by the blast output file.
    :param xml_file: path to the blast output file which should be in xml format
    :return: information provided by blast output file
    """
    curr_info = HSPhit()
    curr_label = ""
    ret = {}
    with open(xml_file, 'r') as file:
        for line in file:
            tag = first_tag_in_string(line)

            if tag == "Iteration_query-def":
                next_hsp = tag_contents(line, tag)
                end = next_hsp.find("_unambig")
                if end != -1:
                    next_hsp = next_hsp[0:end]
                curr_label = next_hsp

            elif tag == "Hsp":
                curr_info = HSPhit()

            elif tag == "/Hsp":
                # Replace element in ret if:
                #       curr_label not in ret
                # OR    current bit-score > ret[*].bit-score
                # OR    bit-score equal
                #   AND     current evalue < ret[*].evalue
                # do NOT replace if bit-score/evalue equal
                replace = not(curr_label in ret
                              and not(curr_info.bit_score > ret[curr_label].bit_score
                                      or(curr_info.bit_score == ret[curr_label].bit_score
                                         and curr_info.evalue < ret[curr_label].evalue)))

                if replace:
                    ret[curr_label] = curr_info

            elif tag == "Hsp_bit-score":
                curr_info.bit_score = float(tag_contents(line, tag))

            elif tag == "Hsp_evalue":
                curr_info.evalue = float(tag_contents(line, tag))

            elif tag == "Hsp_qseq":
                curr_info.qSeq = tag_contents(line, tag)

            elif tag == "Hsp_hseq":
                curr_info.hSeq = tag_contents(line, tag)

            elif tag == "Hsp_midline":
                curr_info.midLine = tag_contents(line, tag)

            elif tag == "Hsp_hit-from":
                curr_info.hitFrom = tag_contents(line, tag)

            elif tag == "Hsp_hit-to":
                curr_info.hitTo = tag_contents(line, tag)

            elif tag == "Hsp_query-from":
                curr_info.queryFrom = tag_contents(line, tag)

            elif tag == "Hsp_query-to":
                curr_info.queryTo = tag_contents(line, tag)

            elif tag == "Hit_def":
                curr_info.hit = tag_contents(line, tag)

    return ret


def first_tag_in_string(line):
    """
    Helper function for the extract alignments function
    :param line: a line in the blast output file
    :return: the information in the line without the tags
    """
    start = line.find("<")
    end = line.find(">", start + 1)

    if end == -1 or start == -1:
        return ""

    return line[start+1:end]


def tag_contents(line, tag):
    """
    Helper function for the extract alignments function
    :param line: a line in the blast output file
    :param tag: a given tag. Information gets extracted after this tag occurs
    :return: the information in the line without the tags
    """
    start = line.find("<" + tag + ">")

    if start == -1:
        return ""

    start += len(tag) + 2
    end = line.find("</" + tag + ">", start)

    if end == -1:
        return ""

    return line[start:end]


# ===========================================================
# Misc.

def get_mtimes(*args):
    """
    Map last-modified time to extension (without leading '.')
    time format: mm/dd/yy hh:mm:ss
    :return: dict with 'ext: mtime' elements
    """
    ret = {}

    for file in args:
        _, ext = os.path.splitext(file)
        ret[ext[1:]] = time.strftime(
            '%m/%d/%y %H:%M:%S', time.localtime(os.path.getmtime(file)))

    return ret


def define_column_positions(header):
    """
    Defines which column in the header of the report file is at which position.
    :param header: header row of the SeqScreen output (report file)
    :return: the positions of the different columns
    """
    if header.startswith("#"):
        header = header[1:]

    global funsocs_end_idx, funsocs_start_idx, funsocs_col_names
    head_split = header.rstrip().split("\t")
    ret = {}
    metadata_cols = {
        "query",
        "taxid",
        "tax_id",
        "go",
        "multi_taxids_confidence",
        "go_id_confidence",
        "size",
        "organism",
        "gene_name",
        "uniprot",
        "bsat_hit",
        "vfdb_hit",
        "flag",
    }
    for i, item in enumerate(head_split):
        ret[item] = i

    if "disable_organ" in ret and "virulence_regulator" in ret:
        start = min(ret["disable_organ"], ret["virulence_regulator"])
        end = max(ret["disable_organ"], ret["virulence_regulator"])
        funsoc_indexes = list(range(start, end + 1))
    else:
        start_after = ret.get("go_id_confidence", ret.get("go", -1)) + 1
        end_before = ret.get("size", len(head_split))
        funsoc_indexes = [
            i for i in range(start_after, end_before)
            if head_split[i] not in metadata_cols
        ]

    if funsoc_indexes:
        funsocs_start_idx = funsoc_indexes[0]
        funsocs_end_idx = funsoc_indexes[-1]
        for idx, index in enumerate(funsoc_indexes):
            funsocs_col_names[idx] = head_split[index]

    return ret


def getAllGoTerms(goTerms):
    if(goTerms):
        result = []
        for go in goTerms:
            result += getParents(go)
    return goTerms+list(set(result)-set(goTerms))
   


def getParents(child):
    global goTerms
    results = []
    results.append(child)
    hold = []
    hold.append(child)
    while len(hold) > 0:
        numOfParents = len(hold)
        for i in range(0, len(hold)):
            if hold[i] in goTerms:
                results += goTerms[hold[i]]
                hold += goTerms[hold[i]]
        hold = hold[numOfParents:]
    return list(set(results))


def optional_existing_path(path, expected_kind, label):
    """
    Resolve an optional report asset, returning None when it is unavailable.
    """
    resolved = os.path.abspath(path)
    exists = os.path.isdir(resolved) if expected_kind == "dir" else os.path.isfile(resolved)
    if exists:
        return resolved
    print(f"Warning: optional report asset not found: {label} ({resolved}). Using fallback output.", file=sys.stderr)
    return None


def main():
    """
    Main function.
    Generates html output which shows threat potential for each sequence.
    Also shows all information about each sequence.
    """

    # ===========================================================
    # define the required arguments
    parser = argparse.ArgumentParser(description='HTML report generation')
    parser.add_argument(
        "-r", "--report", help='Seqscreen report file (.txt).', required=True)
    parser.add_argument(
        "-f", "--fasta", help='Sequence file in fasta format.', required=True)
    parser.add_argument("-b", "--blast_output",
                        help='Blast output xml file.', required=True)
    parser.add_argument(
        "-o", "--output", help='Where to save output archive.', default="result.zip")
    parser.add_argument("-d", "--dependencies",
                        help='Location of folder containing CSS and JS dependencies.', required=True)
    parser.add_argument("-t", "--template",
                        help='Location of template html file.', required=True)
    parser.add_argument("-g", "--go_names",
                        help='Location of go_names tsv file.', required=True)
    parser.add_argument("-n", "--go_network",
                        help='Location of go_network tsv file.', required=True)
    parser.add_argument(
        "--go_template", help='Location of go term visualization html template file.', required=True)
    parser.add_argument("--funsocs",
                        help='Location of funsocs txt file.', required=True)
    parser.add_argument(
        "--version", type=str, help='Version number of seqscreen', required=True)
    parser.add_argument( "--mode", type=str, help='Default/Sensitive Mode', required=True)
    parser.add_argument( "--rflag", type=str, help='Information about the runtime flags', required=True)
    args = parser.parse_args()

    # ===========================================================
    # Process program arguments

    # file format check
    if not args.report.endswith('.tsv'):
       sys.exit("Report file is in invalid format!")
    if not args.blast_output.endswith('.xml'): 
       sys.exit("Blast Output file is in invalid format!")
    if not args.template.endswith('.html'):
       sys.exit("Template file is in invalid format!")
    if not args.go_names.endswith('.txt'):
       sys.exit("go_names file is in invalid format!")
    if not args.go_network.endswith('.txt'):
       sys.exit("go_network file is in invalid format!")
    if not args.go_template.endswith('.html'):
       sys.exit("go_template file is in invalid format!")
    if not args.funsocs.endswith('.txt'):
       sys.exit("funsocs file is in invalid format!")

    #updating the mode as appropriate 
    global mode, num_cols
    if(args.mode == 'fast'):
        mode = 'Default'
    elif args.mode == 'ont':
        mode = 'ONT'
    else:
        mode = 'Sensitive'
    rflag = 'None'
    if(args.rflag != 'None'):
        rflag = re.sub(r'[,]', r' , ',args.rflag) 
        rflag = re.sub(r'[:]', r' : ',rflag)

    # Get optional report presentation assets. Missing legacy assets should not
    # prevent SeqScreen results from being reported.
    dependencies = optional_existing_path(args.dependencies, "dir", "HTML dependencies")
    template_file = optional_existing_path(args.template, "file", "technical HTML template")
    go_template = optional_existing_path(args.go_template, "file", "GO HTML template")
    go_names_file = optional_existing_path(args.go_names, "file", "GO term names")
    go_network_file = optional_existing_path(args.go_network, "file", "GO term network")
    funsocs_file = optional_existing_path(args.funsocs, "file", "functional category descriptions")

    # Get location of output directory
    if args.output:
        archive_name = os.path.abspath(args.output)
    os.makedirs(archive_name, exist_ok=True)
    # Specify output files
    index_file = os.path.join(archive_name, "index.html")
    out_file = os.path.join(archive_name, "out.html")
    # Specify go visualization output file
    go_out_file = os.path.join(archive_name, "go.html")
    # ===========================================================

    # ===========================================================
    # Process inputs

    # extracts important information from the blast_output file
    alignments = extract_alignments(args.blast_output)

    # extracts the sequence names and its lengths
    seq_names_lengths = load_all(args.fasta)

    # parse go_names tsv file into unique dict
    go_names = load_go_names(go_names_file)

    # parse go_network tsv file into dict
    go_network = load_go_network(go_network_file)

    # combines everything to get a list with all sequences and their specific information
    seqs = load_sequences(args.report, seq_names_lengths, alignments, go_names)

    # counts how many seqs are in each funsocs category
    tally = count_funsocs(seqs)

    # Creates a list for the funsocs count table
    funsocs_count = build_funsocs_count(tally)
    
    # Adding description to the funsocs
    funsocs_title = set_title(funsocs_file)
    
    #Create a list of the bpoc names that will become the column names
    funsocs_names = build_funsocs_names(funsocs_count,funsocs_title)

    # format all collected information about the sequences
    # so that they can be assigned to a variable in JavaScript
    formatted_seqs = json.dumps(format_sequence_information(seqs))
    formatted_go_seqs = json.dumps(format_go_information(seqs))

    # get time of last modification of input files
    mtimes = get_mtimes(args.report, args.fasta, args.blast_output)

    # ===========================================================

    # ===========================================================
    # Save outputs

    # copy dependencies
    if dependencies:
        copy_tree(dependencies, archive_name)

    # create the output file with an html template and the calculated variables
    create_executive_report(index_file, seqs, funsocs_count, funsocs_names, args.version, mode, rflag, args.fasta, mtimes)
    if template_file and Environment:
        create_output_file(template_file, out_file, "out",funsocs_count,args.version,funsocs_names)
    else:
        create_fallback_output_file(out_file, funsocs_count, args.version, funsocs_names)

    if go_template and Environment:
        create_go_file(go_template, go_out_file, formatted_go_seqs, go_network)
    else:
        create_fallback_go_file(go_out_file, go_network)
    write_file(seqs,archive_name,mode,rflag,args.fasta,len(seqs),mtimes)
    # ===========================================================

    # removes the created output file after packing it with the other files into a zip archive
    # os.remove(outFile)


if __name__ == "__main__":
    main()
