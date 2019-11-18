#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import pytest
import json
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from basicfunctions import GeneralFunctions as G

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

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")

start = (
        "Create new config files or start pipeline with previously "
        "generated files?\ntype (n)ew or (s)tart:\n")
species = (
        "Please specify target species (comma separated) "
        "or type help for format examples: \n")
path = (
        "Please specify a path to use as the working directory "
        "or hit return to use the current working directory:\n")
targets = (
        "Search for config files for (a)ll or (s)elect targets:\n")

# full_run paths
config_dir = os.path.join(
    "/", "primerdesign", "test", "Lactobacillus_curvatus", "config")
conffile = os.path.join(testfiles_dir, "fullrun_config.json")
testconfig = os.path.join(config_dir, "config.json")
summ_dir = os.path.join(
        "/", "primerdesign", "test", "Summary", "Lactobacillus_curvatus")
genomic_dir = os.path.join(
        "/", "primerdesign", "test", "Lactobacillus_curvatus", "genomic_fna")

def prepare_testfiles():
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
                            record.id +
                            " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbpath

    def create_customblastdb():
        cmd = [
            "makeblastdb", "-in", dbpath, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    G.create_directory(os.path.dirname(dbpath))
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


def prepare_files():
    if os.path.isdir(genomic_dir):
        shutil.rmtree(genomic_dir)
    G.create_directory(genomic_dir)
    G.create_directory(config_dir)
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
        newfilename = ".".join(
                filename.split(".ffn")[0].split("v")) + "_genomic.fna"
        outpath = os.path.join(genomic_dir, newfilename)
        mockfna = SeqRecord(
            Seq("".join(sequences)),
            id="MOCK_" + filename.split("_")[1],
            name="MOCK_" + filename.split("_")[1],
            description="Lactobacillus curvatus")
        with open(outpath, "w") as o:
            SeqIO.write(mockfna, o, "fasta")
    maxcontig = os.path.join(genomic_dir, "GCF_007MOCKv1_genomic.fna")
    with open(maxcontig, "w") as f:
        f.write(">GCF_007MOCKv1_file\n")
        for i in range(0, 500):
            f.write(">GCF_007MOCKv1_" + str(i) + "\n")


def assert_ref_files(nolist=False):
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

    if nolist:
        ref_files.append("potential_specieslist.txt")

    files.sort()
    ref_files.sort()
    assert files == ref_files

def test_run(monkeypatch):
    prepare_testfiles()
    prepare_files()
    os.chdir(os.path.join("/", "primerdesign"))
    from speciesprimer import main
    main(mode="auto")
    maxcontig = os.path.join(genomic_dir, "GCF_007MOCKv1_genomic.fna")
    # compare primer results
    assert os.path.isfile(maxcontig) is False
    assert_ref_files()
    # will not run again
    monkeypatch.setattr('builtins.input', start_oneinput)
    from speciesprimer import main
    main()

    if os.path.isdir(summ_dir):
        shutil.rmtree(summ_dir)

def test_rerun(monkeypatch):
    pandir = os.path.join(
        "/", "primerdesign", "test", "Lactobacillus_curvatus", "Pangenome")
    with open(testconfig) as f:
        for line in f:
            confdict = json.loads(line)
    confdict.update({"qc_gene": ["tuf"]})
    confdict.update({"ignore_qc": True})
    with open(testconfig, "w") as f:
        f.write(json.dumps(confdict))
    if os.path.isdir(pandir):
        shutil.rmtree(pandir)
    monkeypatch.setattr('builtins.input', start_oneinput)
    from speciesprimer import main
    main()
    mfold_dir = os.path.join(pandir, "results", "primer", "mfold")
    mfoldfiles = os.listdir(mfold_dir)
    mfoldfiles.sort()
    assert mfoldfiles == [
            "mfold_failed.csv", "mfold_passed.csv"]
    assert_ref_files()
    if os.path.isdir(summ_dir):
        shutil.rmtree(summ_dir)

def test_repeat_primerQC(monkeypatch):
    with open(testconfig) as f:
        for line in f:
            confdict = json.loads(line)
    confdict.update({"nolist": True})
    with open(testconfig, "w") as f:
        f.write(json.dumps(confdict))
    monkeypatch.setattr('builtins.input', start_oneinput)
    from speciesprimer import main
    main()
    assert_ref_files(nolist=True)


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
        assert os.path.isdir(test) is False

    remove_test_files()

if __name__ == "__main__":
    print(msg)
