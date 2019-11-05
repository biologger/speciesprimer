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
import csv
from Bio import SeqIO

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
#sys.path.append(BASE_PATH.split("tests")[0])
sys.path.append(pipe_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G

testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")

confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "primerdesign", "test"),
    "probe": False, "exception": None, "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["rRNA"],
    "blastdbv5": False, "intermediate": True, "nontargetlist": ["Lactobacillus sakei"]}

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

@pytest.fixture
def config():
    from speciesprimer import CLIconf
    args = AttrDict(confargs)
    nontargetlist = H.create_non_target_list(args.target)
    config = CLIconf(
            args.minsize, args.maxsize, args.mpprimer, args.exception,
            args.target, args.path, args.intermediate,
            args.qc_gene, args.mfold, args.skip_download,
            args.assemblylevel, nontargetlist,
            args.skip_tree, args.nolist, args.offline,
            args.ignore_qc, args.mfethreshold, args.customdb,
            args.blastseqs, args.probe, args.blastdbv5)

    config.save_config()

    return config

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
        skip_download: "n",
        assemblylevel: "complete",
        customdb: "",
        maxseqs: "1000",
        qc_genes: "tuf",
        exception: 'Bacterium_curvatum',
        minsize: "70",
        maxsize: "200",
        designprobe: "y",
        mfold_th: "",
        mpprimer_th: "",
        mfeprimer_th: "",
        ignore_qc: "n",
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

def test_batchassist(config, monkeypatch):   
    from speciesprimer import Config
    from speciesprimer import CLIconf
    
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
    
        assert config.assemblylevel == ["complete"]
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
    
        assert config.assemblylevel == ["complete"]
        assert config.blastseqs == 1000
        assert config.exception == 'Bacterium_curvatum'
        assert config.maxsize == 200
    
    
    
    