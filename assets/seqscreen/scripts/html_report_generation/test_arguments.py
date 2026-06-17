#!/usr/bin/env python3
"""
Argument tests for html_report_generation
"""
import os
import subprocess
import pytest

# when -o is not specified, result.zip will be placed in current directory
DEFAULT_OUTPUT_PATH = "./result.zip"
# Working Directory == directory where this file is located
WD = os.path.dirname(os.path.abspath(__file__))
# all the argument's paths
REPORT = WD + "/data/seqscreen_report.tsv"
FASTA = WD + "/data/seqs.fasta"
BLASTX = WD + "/data/seqs.fasta.ur100.xml"
OUT = WD + "/result.zip"
LIBS = WD + "/libs/"
TEMPLATE = WD + "/data/template.html"
GO_NAMES = WD + "/data/go_names.txt"
GO_NETWORK = WD + "/data/go_network.txt"
GO_TEMPLATE = WD + "/data/go_template.html"


INVALID_INPUTS = [
    # missing required argument
    ((f"-f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} --go_template {GO_TEMPLATE}"), 2),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK}"), 2),
    # invalid file endings: report, fasta, blast_output, template, go_names, go_network, go_template
    ((f"-r {os.path.splitext(REPORT)[0] + '.csv'} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {os.path.splitext(FASTA)[0] + '.xml'} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {FASTA} -b {os.path.splitext(BLASTX)[0] + '.json'} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {os.path.splitext(TEMPLATE)[0] + '.txt'} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {os.path.splitext(GO_NAMES)[0] + '.csv'} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {os.path.splitext(GO_NETWORK)[0] + '.csv'} --go_template {GO_TEMPLATE}"), 1),
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {os.path.splitext(GO_TEMPLATE)[0] + '.txt'}"), 1)
]
VALID_INPUTS = [
    # all arguments
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -o {OUT} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 0),
    # all required arguments
    ((f"-r {REPORT} -f {FASTA} -b {BLASTX} -d {LIBS} -t {TEMPLATE} -g {GO_NAMES} -n {GO_NETWORK} --go_template {GO_TEMPLATE}"), 0)
]

START_SCRIPT = WD + "/generateHtmlReport.py"


# ===========================================================
# Argument tests
@pytest.mark.parametrize("params, expected_status", INVALID_INPUTS + VALID_INPUTS)
def test_return_status(params, expected_status):
    """
    Test the return status for invalid and valid inputs
    """
    status = subprocess.call(f'python3 {START_SCRIPT} {params}', shell=True)
    assert expected_status == status

# Output file tests
@pytest.mark.parametrize("params, expected_status", INVALID_INPUTS + VALID_INPUTS)
def test_output(params, expected_status):
    """
    Test if a output folder is created
    Should create folder for valid inputs, should not create folder for invalid inputs
    """
    # output path is either DEFAULT_OUTPUT_PATH or whatever is specified by '-o'
    output_path = ""
    if "-o " in params:
        split_params = params.split(' ')
        output_path = split_params[split_params.index('-o') + 1]
    else:
        output_path = DEFAULT_OUTPUT_PATH


    if os.path.isfile(output_path):
        os.remove(output_path)

    subprocess.call(f'python3 {START_SCRIPT} {params}', shell=True)

    if expected_status == 0:
        assert os.path.isfile(output_path)
    else:
        assert not os.path.isfile(output_path)
