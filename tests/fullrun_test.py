#!/usr/bin/env python3
# -*- coding: utf-8 -*-

msg = (
"""
Works only in the Docker container!!!
- Start the container
    sudo docker start {Containername}
- Start an interactive terminal in the container
    sudo docker exec -it {Containername} bash
- Start the test in the container terminal
    pytest -vv /pipeline/tests/fullrun_test.py
"""
)

import os
import sys
import shutil
import pytest
import json
import time
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")

#            import traceback
#            traceback.print_exc()

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G

testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")

start = ("Create new config files or start pipeline with previously "
        "generated files?\ntype (n)ew or (s)tart:\n")
species = ("Please specify target species (comma separated) "
        "or type help for format examples: \n")
path = ("Please specify a path to use as the working directory "
            "or hit return to use the current working directory:\n")
targets = (
    "Search for config files for (a)ll or (s)elect targets:\n")

def prepare_testfiles():
    G.create_directory(os.path.dirname(dbpath))
    def prepare_tmp_db():
        t = os.path.join(BASE_PATH, "testfiles", "tmp_config.json")
        # Docker only
        tmp_path = os.path.join(pipe_dir, "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        shutil.copy(t, tmp_path)

    def change_tmp_db():
        tmp_path = os.path.join(pipe_dir, "tmp_config.json")
        with open(tmp_path) as f:
            for line in f:
                tmp_dict = json.loads(line)
        tmp_dict["new_run"].update({'modus': "continue", "targets": None})
        with open(tmp_path, "w") as f:
            f.write(json.dumps(tmp_dict))

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        with open(dbpath, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbpath

    def create_customblastdb():
        cmd = [
            "makeblastdb", "-in", dbpath, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    dbinputfiles()
    create_customblastdb()
    prepare_tmp_db()
    change_tmp_db()

def start_oneinput(prompt):
    prompt_dict = {
        start: "s",
        targets: "s",
        species: "Lactobacillus curvatus",
        path: "/primerdesign/test"
    }
    val = prompt_dict[prompt]
    return val

def test_run(monkeypatch):
    def prepare_files():
        genomic_dir = os.path.join(
                "/", "primerdesign", "test", "Lactobacillus_curvatus", "genomic_fna")
        if os.path.isdir(genomic_dir):
            shutil.rmtree(genomic_dir)
        G.create_directory(genomic_dir)
        config_dir = os.path.join(
            "/", "primerdesign", "test", "Lactobacillus_curvatus", "config")
        G.create_directory(config_dir)
        conffile = os.path.join(testfiles_dir, "fullrun_config.json")
        testconfig = os.path.join(config_dir, "config.json")
        shutil.copy(conffile, testconfig)
        files = ["GCF_001981925v1_20190923.ffn", "GCF_003410375v1_20190923.ffn"]
        ffn_files_dir = os.path.join(testfiles_dir, "ffn_files")

        for filename in files:
            filepath = os.path.join(ffn_files_dir, filename)
            sequences = []
            with open(filepath) as f:
                records = list(SeqIO.parse(f, "fasta"))
            for record in records:
                seq = str(record.seq)
                sequences.append(seq)
            newfilename = ".".join(filename.split(".ffn")[0].split("v")) + "_genomic.fna"
            outpath = os.path.join(genomic_dir, newfilename)
            mockfna = SeqRecord(
                Seq("".join(sequences)),
                id="MOCK_" + filename.split("_")[1],
                name="MOCK_" + filename.split("_")[1],
                description="Lactobacillus curvatus")
            with open(outpath, "w") as o:
                SeqIO.write(mockfna, o, "fasta")

    prepare_testfiles()
    prepare_files()
    os.chdir(os.path.join("/", "primerdesign"))
    monkeypatch.setattr('builtins.input', start_oneinput)
    from speciesprimer import main
    main(mode="auto")
    summ_dir = os.path.join(
                "/", "primerdesign", "test", "Summary", "Lactobacillus_curvatus")
    files = []
    for filename in os.listdir(summ_dir):
        files.append(filename)

    today = time.strftime("%Y_%m_%d", time.localtime())
    configfile = "config_" + today + ".json"
    stats = "Lb_curva_pipeline_stats_" + today + ".txt"

    ref_files = [
        "core_gene_alignment.aln", "Lb_curva_tree.nwk",
        "Lb_curva_primer.csv", "Lb_curva_qc_sequences_details.csv",
        "mostcommonhits.csv", "Lb_curva_qc_sequences.csv", configfile,
        stats]

    files.sort()
    ref_files.sort()
    assert files == ref_files

    # open config file and change keep intermedaite files to False
    # rerun on exisitng files
    # change also probe option
    # add excluded gis
    # test parallel functions non-parallel to get more covered lines
    # add no list option
    # include argparser/commandline function?
    # include errors to provoke an error report (multiple targets)

def test_end():
    def remove_test_files():
        test =  os.path.join("/", "primerdesign", "test")
        shutil.rmtree(test)
        tmp_path = os.path.join("/", "pipeline", "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.chdir(BASE_PATH)
        assert os.path.isdir(test) == False

    remove_test_files()

if __name__ == "__main__":
    print(msg)
