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
    pytest -vv /pipeline/tests/speciesprimer_test.py
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


# prompts
start = ("Create new config files or start pipeline with previously "
        "generated files?\ntype (n)ew or (s)tart:\n")
species = ("Please specify target species (comma separated) "
        "or type help for format examples: \n")
path = ("Please specify a path to use as the working directory "
            "or hit return to use the current working directory:\n")
skip_tree = (
    "Skip the core gene alignment and tree for visualization and "
    "troubleshooting.\n(y)es/(n)o, default=(n)\n> ")
offline = (
    "Work offline with local genome assemblies?"
    "\n(y)es/(n)o, default=(n)\n> ")
skip_download = (
    "Skip the download of Genomes from NCBI?"
    "\n(y)es/(n)o, default=(n)\n> ")
assemblylevel = (
    "Limit downloads of Genomes to assembly status\n"
    "options: complete, chromosome, scaffold, contig, "
    "all (comma separated), default=all\n> ")
customdb = (
    "Do you want to use a custom database for blastn?\n"
    "Specifiy the absolute filepath of the custom database "
    "(e.g. '/home/blastdb/nontarget.fasta') or hit return"
    " to skip, default=None\n> ")
maxseqs = (
    "Set the number of sequences per BLAST search. "
    "Decrease the number of sequences if BLAST slows down "
    "due to low memory."
    "\noptions = [100, 500, 1000, 2000, 5000]"
    ", default=1000\n> ")
qc_genes = (
    "Gene(s) (comma separated) for BLAST in the initial quality "
    "control step.\noptions: rRNA, tuf, recA, dnaK, pheS"
    ", default=rRNA\n> ")
exception = (
    "Primer binding to this non-target species is tolerated.\n"
    "Provide a species name or hit return to skip:\n> ")
minsize = ("Minimal Amplicon size\ndefault=70\n> ")
maxsize = ("Maximal amplicon size\ndefault=200\n> ")
designprobe = (
    "Do you want primer3 to design an internal probe?"
    "[Experimental!]"
    "\n(y)es/(n)o, default=(n)\n> ")
mfold_th = (
    "Delta G threshold for secondary structures in PCR products"
    " at 60 degree Celsius calculated by mfold\ndefault=-3.0\n> ")
mpprimer_th = (
    "Mpprimer threshold for delta G values calculated by "
    "MPprimer_dimer_check.pl\ndefault=-3.5\n> ")
mfeprimer_th = (
    "MFEprimer threshold for nontarget sequence PPC "
    "higher values mean more stringent selection.\ndefault=90\n> ")
ignore_qc = (
    "Do you want to include genomes that did not"
    " pass quality control?\ndefault=(n)\n> ")
blastdbv5 = (
    "Do you have the Version 5 of the BLAST DB? \ndefault=(n)\n> ")
intermediate = (
    "Do you want to keep intermediate files?\ndefault=(n)\n> ")
nolist = (
    "Do you want to perform specificity check without the "
    "(non-target) species list (for all sequences in the DB)?"
    "\nNot recommended for nt DB! May be used with a custom DB"
    "\ndefault=(n)\n> ")
forall = (
    "Use this value for all targets?\n(y)es/(n)o, default=(y)\n> ")

targets = (
    "Search for config files for (a)ll or (s)elect targets:\n")

def good_input(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus curvatus",
        path: "/primerdesign/test",
        skip_tree: "n",
        offline: "n",
        skip_download: "y",
        assemblylevel: "complete",
        customdb: dbpath,
        maxseqs: "1000",
        qc_genes: "tuf",
        exception: 'Bacterium_curvatum',
        minsize: "70",
        maxsize: "200",
        designprobe: "y",
        mfold_th: "",
        mpprimer_th: "",
        mfeprimer_th: "",
        ignore_qc: "y",
        blastdbv5: "n",
        intermediate: "y",
        nolist: "n",
        forall: "y"
    }
    val = prompt_dict[prompt]
    return val

def start_input(prompt):
    prompt_dict = {
        start: "s",
        targets: "s",
        species: "Lactobacillus curvatus",
        path: "/primerdesign/test"
    }
    val = prompt_dict[prompt]
    return val

def bad_input(prompt):
    prompt_dict = {
        start: "n",
        species: "",
        path: "",
        skip_tree: "",
        offline: "",
        skip_download: "",
        assemblylevel: "",
        customdb: "",
        maxseqs: "",
        qc_genes: "",
        exception: "",
        minsize: "",
        maxsize: "",
        designprobe: "",
        mfold_th: "",
        mpprimer_th: "",
        mfeprimer_th: "",
        ignore_qc: "",
        blastdbv5: "",
        intermediate: "",
        nolist: "",
        forall: ""
    }
    val = prompt_dict[prompt]
    return val

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
            cmd, printcmd=False, logcmd=False, log=False, printoption=False)

    dbinputfiles()
    create_customblastdb()
    prepare_tmp_db()
    change_tmp_db()


def test_batchassist(monkeypatch):
    from speciesprimer import Config
    from speciesprimer import CLIconf

    test =  os.path.join("/", "primerdesign", "test")
    if os.path.isdir(test):
        shutil.rmtree(test)
    prepare_testfiles()
    monkeypatch.setattr('builtins.input', good_input)
    conf_from_file = Config()
    targets = conf_from_file.get_targets()
    nontargetlist = []
    for target in targets:
        (
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist,
            offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5
                ) = conf_from_file.get_config(target)

        config = CLIconf(
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, nontargetlist, skip_tree,
            nolist, offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5)

        assert config.assemblylevel == ["offline"]
        assert config.blastseqs == 1000
        assert config.exception == 'Bacterium_curvatum'
        assert config.maxsize == 200

        config.save_config()

    monkeypatch.setattr('builtins.input', start_input)
    conf_from_file = Config()
    for target in targets:
        (
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist,
            offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5
                ) = conf_from_file.get_config(target)

        config = CLIconf(
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, nontargetlist, skip_tree,
            nolist, offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5)

        assert config.assemblylevel == ["offline"]
        assert config.blastseqs == 1000
        assert config.exception == 'Bacterium_curvatum'
        assert config.maxsize == 200

def test_run(monkeypatch):

    def prepare_files():
        genomic_dir = os.path.join(
                "primerdesign", "test", "Lactobacillus_curvatus", "genomic_fna")
        if os.path.isdir(genomic_dir):
            shutil.rmtree(genomic_dir)
        G.create_directory(genomic_dir)
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

    prepare_files()
    os.chdir(os.path.join("/", "primerdesign"))

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
