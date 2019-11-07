#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from speciesprimer import Config
from speciesprimer import CLIconf

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
qc_gene = (
    "Gene(s) (comma separated) for BLAST in the initial quality "
    "control step.\noptions: rRNA, tuf, recA, dnaK, pheS"
    ", default=rRNA\n> ")
exception = (
    "Primer binding to this non-target species is tolerated.\n"
    "Provide a species name or hit return to skip:\n> ")
minsize = ("Minimal Amplicon size\ndefault=70\n> ")
maxsize = ("Maximal amplicon size\ndefault=200\n> ")
probe = (
        "Do you want primer3 to design an internal probe? [Experimental!]"
        "\n(y)es/(n)o, default=(n)\n> ")
mfold = (
    "Delta G threshold for secondary structures in PCR products"
    " at 60 degree Celsius calculated by mfold\ndefault=-3.0\n> ")
mpprimer = (
    "Mpprimer threshold for delta G values calculated by "
    "MPprimer_dimer_check.pl\ndefault=-3.5\n> ")
mfethreshold = (
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

def alldef_input(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus",
        path: "",
        skip_tree: "",
        offline: "",
        skip_download: "",
        assemblylevel: "",
        customdb: "",
        maxseqs: "",
        qc_gene: "",
        exception: '',
        minsize: "",
        maxsize: "",
        probe: "",
        mfold: "",
        mpprimer: "",
        mfethreshold: "",
        ignore_qc: "",
        blastdbv5: "",
        intermediate: "",
        nolist: "",
        forall: ""
    }
    val = prompt_dict[prompt]
    return val

def nodef_input():
    pass

def good_input(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus helveticus",
        path: "/primerdesign/test",
        skip_tree: "n",
        offline: "n",
        skip_download: "y",
        assemblylevel: "complete",
        customdb: dbpath,
        maxseqs: "1000",
        qc_gene: "tuf",
        exception: 'Bacterium_curvatum',
        minsize: "70",
        maxsize: "200",
        probe: "y",
        mfold: "",
        mpprimer: "",
        mfethreshold: "",
        ignore_qc: "y",
        blastdbv5: "n",
        intermediate: "y",
        nolist: "n",
        forall: "n"
    }
    val = prompt_dict[prompt]
    return val

def start_input(prompt):
    prompt_dict = {
        start: "s",
        targets: "a",
        path: "/primerdesign/test"
    }
    val = prompt_dict[prompt]
    return val

def fail_startinput(prompt):
    prompt_dict = {
        start: "s",
        targets: "q",
        path: "/primerdesign/test"
    }
    val = prompt_dict[prompt]
    return val

def start_wronginput(prompt):
    prompt_dict = {
        start: "s",
        targets: "s",
        species: "Lactobacillus bifermentans",
        path: "/primerdesign/test"
    }
    val = prompt_dict[prompt]
    return val

def start_oneinput(prompt):
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
        skip_tree: "q",
        offline: "q",
        skip_download: "q",
        assemblylevel: "alle",
        customdb: "customdb.fas",
        maxseqs: "66",
        qc_gene: "tufrrna",
        exception: "Lactobacillus_curvatus, Lactobacillus helveticus",
        minsize: "50.36",
        maxsize: "twohundred",
        probe: "q",
        mfold: "minusthree",
        mpprimer: "minusthreepointfive",
        mfethreshold: "99.9",
        ignore_qc: "q",
        blastdbv5: "q",
        intermediate: "q",
        nolist: "q",
        forall: "q"
    }
    val = prompt_dict[prompt]
    return val

def bad_input2(prompt):
    prompt_dict = {
        start: "n",
        species: "help",
        path: "",
        skip_tree: "q",
        offline: "q",
        skip_download: "q",
        assemblylevel: "alle",
        customdb: "customdb.fas",
        maxseqs: "66",
        qc_gene: "tufrrna",
        exception: "Lactobacillus_curvatus, Lactobacillus helveticus",
        minsize: "50.36",
        maxsize: "twohundred",
        probe: "q",
        mfold: "minusthree",
        mpprimer: "minusthreepointfive",
        mfethreshold: "99.9",
        ignore_qc: "q",
        blastdbv5: "q",
        intermediate: "q",
        nolist: "q",
        forall: "q"
    }
    val = prompt_dict[prompt]
    return val

def fail_input(prompt):
    if prompt:
        return "q"

class batchassist_mock():
    def __init__(self, good_input, bad_input):
        self.repeat = 0
        self.good_input = good_input
        self.bad_input = bad_input
    def prompt_input(self, prompt):
        if self.repeat == 0:
            result = self.bad_input(prompt)
            self.repeat +=1
        else:
            result = self.good_input(prompt)
            self.repeat = 0
        return result

def get_config_from_file(conf_from_file):
    configfilepaths = []
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

        config.save_config()

        cfilepath = os.path.join(config.path, config.target, "config", "config.json")
        configfilepaths.append(cfilepath)
    return configfilepaths


def prepare_testfiles():
    G.create_directory(os.path.dirname(dbpath))
    def prepare_tmp_db():
        t = os.path.join(testfiles_dir, "tmp_config.json")
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

def test_sys_exit(monkeypatch):
    monkeypatch.setattr('builtins.input', fail_input)
    with pytest.raises(SystemExit):
        conf_from_file = Config()

    monkeypatch.setattr('builtins.input', fail_startinput)
    with pytest.raises(SystemExit):
        conf_from_file = Config()


# mock >primer_csv = os.path.join(path, "Summary", target, abbr + "_primer.csv")

def test_batchassist(monkeypatch, caplog):


    caplog.set_level(logging.INFO)
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


    monkeypatch.setattr('builtins.input', start_input)
    conf_from_file = Config()
    targets = conf_from_file.get_targets()
    reftarget = ["Lactobacillus_curvatus", "Lactobacillus_helveticus"]
    targets.sort()
    reftarget.sort()
    for i, target in enumerate(targets):
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

        assert config.target == reftarget[i]
        assert config.assemblylevel == ["offline"]
        assert config.blastseqs == 1000
        assert config.exception == 'Bacterium_curvatum'
        assert config.maxsize == 200


    monkeypatch.setattr('builtins.input', start_oneinput)
    conf_from_file = Config()
    targets = conf_from_file.get_targets()

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

        assert config.target == "Lactobacillus_curvatus"
        assert config.assemblylevel == ["offline"]
        assert config.blastseqs == 1000
        assert config.exception == 'Bacterium_curvatum'
        assert config.maxsize == 200

def test_default_input(monkeypatch):
    defaultinputref = os.path.join(testfiles_dir, "default_config.json")
    with open(defaultinputref) as f:
        for line in f:
            defaultinputref_dict = json.loads(line)
    monkeypatch.setattr('builtins.input', alldef_input)




def test_wrong_input(monkeypatch):
    mock = batchassist_mock(good_input, bad_input)
    monkeypatch.setattr('builtins.input', mock.prompt_input)
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

        config.save_config()

    mock = batchassist_mock(start_wronginput, start_oneinput)
    monkeypatch.setattr('builtins.input', mock.prompt_input)
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

        config.save_config()
