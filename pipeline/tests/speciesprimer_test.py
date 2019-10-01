#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import shutil
import pytest


BASE_PATH = os.path.abspath(__file__).split("tests")[0]
new_path = os.path.join(BASE_PATH, "bin")
sys.path.append(new_path)

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G

confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "home", "primerdesign", "test"),
    "probe": False, "exception": None, "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["tuf"],
    "blastdbv5": True, "intermediate": False, "nontargetlist": []}

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

def test_CLIconf(config):

    assert config.minsize == confargs['minsize']
    assert config.maxsize == confargs['maxsize']
    assert config.ignore_qc == confargs['ignore_qc']
    assert config.mfethreshold == confargs['mfethreshold']
    assert config.target == confargs['target']
    assert config.nolist == confargs['nolist']
    assert config.skip_tree == confargs['skip_tree']
    assert config.blastseqs == confargs['blastseqs']
    assert config.mfold == confargs['mfold']
    assert config.mpprimer == confargs['mpprimer']
    assert config.offline == confargs['offline']
    assert config.probe == confargs['probe']
    assert config.exception == confargs['exception']
    assert config.customdb == confargs['customdb']
    assert config.skip_download == confargs['skip_download']
    assert config.assemblylevel == confargs['assemblylevel']
    assert config.qc_gene == confargs['qc_gene']
    assert config.blastdbv5 == confargs['blastdbv5']

def test_auto_run_config():
    from speciesprimer import auto_run
    t = os.path.join(BASE_PATH, "tests", "testfiles", "tmp_config.json")
    # Docker only
    tmp_path = os.path.join("/", "home", "pipeline", "tmp_config.json")
    shutil.copy(t, tmp_path)
    targets, conf_from_file, use_configfile = auto_run()

    assert targets == ["Lactobacillus_curvatus"]
    assert use_configfile == True
    for target in targets:
        target = target.capitalize()
        if use_configfile:
            (
                minsize, maxsize, mpprimer, exception, target, path,
                intermediate, qc_gene, mfold, skip_download,
                assemblylevel, skip_tree, nolist,
                offline, ignore_qc, mfethreshold, customdb,
                blastseqs, probe, blastdbv5
            ) = conf_from_file.get_config(target)

    assert minsize == confargs['minsize']
    assert maxsize == confargs['maxsize']
    assert ignore_qc == confargs['ignore_qc']
    assert mfethreshold == confargs['mfethreshold']
    assert target == confargs['target']
    assert nolist == confargs['nolist']
    assert skip_tree == confargs['skip_tree']
    assert blastseqs == confargs['blastseqs']
    assert mfold == confargs['mfold']
    assert mpprimer == confargs['mpprimer']
    assert offline == confargs['offline']
    assert probe == confargs['probe']
    assert exception == confargs['exception']
    assert customdb == confargs['customdb']
    assert skip_download == confargs['skip_download']
    assert assemblylevel == confargs['assemblylevel']
    assert qc_gene == confargs['qc_gene']
    assert blastdbv5 == confargs['blastdbv5']

def clean_before_tests(config):
    shutil.rmtree(os.path.join(config.path, config.target))

def clean_up_after_test():
    tmp_path = os.path.join("/", "home", "pipeline", "tmp_config.json")
    os.remove(tmp_path)

def test_DataCollection(config, printer):
    clean_before_tests(config)
    from speciesprimer import DataCollection
    DC = DataCollection(config)

    def test_get_email_from_config(config):
        email = DC.get_email_for_Entrez()
        assert email == "biologger@protonmail.com"
        return email
    
    def internet_connection(printer):
        from urllib.request import urlopen
        try:
            response = urlopen('https://www.google.com/', timeout=5)
            return True
        except: 
            print("No internet connection!!! Skip online tests")
            return False    

    def test_get_taxid(config):
        email = DC.get_email_for_Entrez()
        taxid, email = DC.get_taxid(config.target)
        assert taxid == '28038'
        clean_up_after_test()
        return taxid
    
    def test_ncbi_download(taxid, email):
        DC.get_ncbi_links(taxid, email, 1)
        DC.ncbi_download()

    def test_prokka_is_installed():
        cmd = "prokka --citation"
        lines = G.read_shelloutput(cmd, printcmd=False, logcmd=False, printoption=False)
        assert lines[0] == "If you use Prokka in your work, please cite:"
        
    def prepare_prokka(config):
        testdir = os.path.join(BASE_PATH, "tests", "testfiles")
        targetdir = os.path.join(config.path, config.target)
        fileformat = ["fna", "gff", "ffn"]
        for fo in fileformat:
            for files in os.listdir(testdir):
                if files.endswith("."+fo):
                    fromfile = os.path.join(testdir, files)
                    tofile =  os.path.join(targetdir, fo + "_files", files)
                    shutil.copy(fromfile, tofile)
                    if fo == "fna":
                        tofile = os.path.join(targetdir, "genomic_fna", files)
                        shutil.copy(fromfile, tofile)
                    
    def test_run_prokka():
        annotation_dirs, annotated = DC.run_prokka()
        assert annotated == ["GCF_004088235v1"]
        DC.copy_genome_files()
        
    email = test_get_email_from_config(config)
    DC.prepare_dirs()
    if internet_connection(printer):
        taxid = test_get_taxid(config)
        test_ncbi_download(taxid, email)
    test_prokka_is_installed()
    prepare_prokka(config)
    test_run_prokka()
    
def test_roary_is_installed():
    cmd = "roary --w"
    lines = G.read_shelloutput(cmd, printcmd=False, logcmd=False, printoption=False)
    assert lines[1] == "3.12.0"
