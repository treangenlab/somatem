# The following libraries are used:

Bootstrap v4.3.1 (https://getbootstrap.com/docs/4.3/getting-started/download/)

jQuery v3.3.1 (https://jquery.com/download/)

# Command-line arguments

* -r --report (required)
    * S2Fast report file (.txt)
* -f --fasta (required)
    * Sequence file in fasta format
* -b --blast_output (required) 
    * Blast output xml file
* -o --output (default: result.zip)
    * Output archive destination
* -d --dependencies (required)
    * Folder containing CSS/JS dependencies
* -t --template (required)
    * Main HTML template file
* -g --go_names (required)
    * go_names tsv file
    * should use scripts/html_report_generation/data/go_names.txt
* -n --go_network (required)
    * go_network tsv file
    * should use scripts/html_report_generation/data/go_network.txt
* --go_template (required)
    * GO Term visualization HTML template file
    * should use scripts_html_report_generation/data/go_template.html

# Example command to run
```bash
python generateHtmlReport.py -r path/to/seqscreen_report.txt -f path/to/M8_motif_sequences.fasta -b path/to/M8_motif_sequences.fasta.ur100.xml -d path/to/libs -t path/to/template.html -o output/path/of/choice.zip
```
