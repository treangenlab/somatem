import importlib
import importlib.util
import os, sys, string, subprocess, re, argparse

def getCommandOutput(theCommand, checkForStderr):
    p = subprocess.Popen(theCommand, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    #print checkStderr
    if checkForStderr and checkStderr != "":
       return ""
    else:
       return checkStdout.strip().decode('utf-8')

def print_rg(msg, okay=True):
    startcolor = '\033[92m' if okay else '\033[91m'
    print(startcolor + msg + '\033[0m')

def main():
    parser = argparse.ArgumentParser(description='SeqScreen Environment Checker')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-d', "--databases", help="Path to SeqScreen databases", required=True)
    args = parser.parse_args()

    error=False
    print("Checking python >= 3.6...", end=' ')
    if (sys.version_info[0] < 3) or (sys.version_info[0] == 3 and sys.version_info[1] < 6):
        error_exist = True
        print_rg("\nERROR\tPython version is %s."%(sys.version), False)
        raise
    else:
        print(u'\u2713')


    for library in ["scipy", "sklearn", "jinja2", "Bio", "bitarray"]:
        print("Checking for {} package...".format(library), end=' ')
        lib_spec = importlib.util.find_spec(library)
        if not lib_spec:
            print_rg("\nERROR\t{} was not found".format(library), False)
            error = True
        else:
            print(u'\u2713')

    tools = ["nextflow", "blastx", "blastn", "diamond", "hmmscan", "bowtie2", "time", "centrifuge"]
    for program in tools:
        print("Checking if %s is installed..."%(program), end=' ')
        which_program = getCommandOutput("which %s"%program, False)
        if which_program == "":
            error_exist = True
            print_rg("\nERROR\t%s is not found or path is not properly set!"%(program), False)
            error = True
        else:
            print(u'\u2713')

    databases = [
        ("bowtie2/blacklist.seqs.nt.1.bt2", "bowtie2"), 
        ("rapsearch2/blacklist.seqs.aa", "rapsearch"), 
        ("hmmscan/Pfam-A.hmm", "Pfam"),
        ("diamond/uniref.mini.dmnd", "DIAMOND"),
        ("diamond_1_2/uniref_1_2.dmnd", "DIAMOND_1_2"),
        ("rebase/rebase.fna", "rebase"), 
        ("megares/megares_full.nhr", "megares"),
        ("blast/nt/nt.00.nhd", "blast nt"),
        ("taxonomy/taxa_lookup.txt", "taxa lookup"),
        ("blast/UNIREF100.mini.00.phr", "blast Uniref100"),
        ("annotation_scores.pck", "Uniprot annotation scores"),
        ("funsocs.pck", "Funsocs"),
        ("reference_inference/taxid2seqid.pickle", "reference inference taxid pickle"),
        ("reference_inference/taxa.sqlite", "reference inference taxa ete3 db"),
        ("funsocs_commensal_list.pck", "Uniprots with no funsocs"),
        ("go/go_network.txt", "go network"),
        ("uniref_multitax.pickle", "Uniprot multitax pickle")]
    for f, db_name in databases:
        print("Checking for {} database...".format(db_name), end=' ')
        if not os.path.isfile(os.path.join(args.databases, f)):
            print_rg("\nERROR\tCannot find {} database in {}".format(db_name, f), False)
            error = True
        else:
            print(u'\u2713')
    if not error:
        print_rg("Success! All requirements installed")
    else:
        print_rg("Error when checking for requirements. Please see output above", False)



if __name__ == "__main__":
    main()
