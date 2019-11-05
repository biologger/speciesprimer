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
    pytest -vv /home/pipeline/tests/speciesprimer_fulltest.py
"""
)


import os
import sys
import shutil
import pytest
import json
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASE_PATH = os.path.abspath(__file__).split("tests")[0]
new_path = os.path.join(BASE_PATH, "bin")
sys.path.append(new_path)

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G

testfiles_dir = os.path.join(BASE_PATH, "tests", "testfiles")
ref_data = os.path.join(BASE_PATH, "tests", "testfiles", "ref")

confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "home", "primerdesign", "test"),
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

def clean_before_tests(config):
    shutil.rmtree(os.path.join(config.path, config.target))

def test_DataCollection(config):
    t = os.path.join(BASE_PATH, "tests", "testfiles", "tmp_config.json")
    # Docker only
    tmp_path = os.path.join("/", "home", "pipeline", "tmp_config.json")
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)
    shutil.copy(t, tmp_path)
    clean_before_tests(config)
    from speciesprimer import DataCollection
    DC = DataCollection(config)

    def test_get_email_from_config(config):
        email = DC.get_email_for_Entrez()
        assert email == "biologger@protonmail.com"
        return email

    def internet_connection():
        from urllib.request import urlopen
        try:
            response = urlopen('https://www.google.com/', timeout=5)
            return True
        except:
            print("No internet connection!!! Skip online tests")
            return False

    def test_get_taxid(config):
        taxid, email = DC.get_taxid(config.target)
        assert taxid == '28038'
        return taxid

    def test_ncbi_download(taxid, email):
        DC.get_ncbi_links(taxid, email, 1)
        DC.ncbi_download()
        # clean up
        fna = os.path.join(config.path, config.target, "genomic_fna")
        shutil.rmtree(fna)
        G.create_directory(fna)

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

    def remove_prokka_testfiles():
        fileformat = ["fna", "gff", "ffn"]
        targetdir = os.path.join(config.path, config.target)
        for form in fileformat:
            dirpath = os.path.join(targetdir, form + "_files")
            if os.path.isdir(dirpath):
                shutil.rmtree(dirpath)


    email = test_get_email_from_config(config)
    DC.prepare_dirs()
    if internet_connection():
        taxid = test_get_taxid(config)
        test_ncbi_download(taxid, email)
    else:
        G.create_directory(DC.gff_dir)
        G.create_directory(DC.ffn_dir)
        G.create_directory(DC.fna_dir)
    test_prokka_is_installed()
    prepare_prokka(config)
    test_run_prokka()
    remove_prokka_testfiles()
    G.create_directory(DC.gff_dir)
    G.create_directory(DC.ffn_dir)
    G.create_directory(DC.fna_dir)

def test_prokka_and_roary(config):
    from speciesprimer import DataCollection
    from speciesprimer import PangenomeAnalysis
    genomic_dir = os.path.join(config.path, config.target, "genomic_fna")

    def remove_prokka_testfiles():
        fileformat = ["fna", "gff", "ffn"]
        targetdir = os.path.join(config.path, config.target)
        for form in fileformat:
            dirpath = os.path.join(targetdir, form + "_files")
            if os.path.isdir(dirpath):
                shutil.rmtree(dirpath)

    def prepare_files():
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
                description=" ".join(config.target.split("_")))
            with open(outpath, "w") as o:
                SeqIO.write(mockfna, o, "fasta")

    prepare_files()
    DC = DataCollection(config)
    PA = PangenomeAnalysis(config)
    DC.collect()
    today = time.strftime("%Y%m%d", time.localtime())
    dirnames = ["GCF_001981925v1_" + today, "GCF_003410375v1_" + today]
    PA.run_pangenome_analysis()

    for dirname in dirnames:
        dirpath = os.path.join(PA.target_dir, dirname)
        filelist = []
        for files in os.listdir(dirpath):
            filelist.append(files)

        assert len(filelist) == 12
        shutil.rmtree(dirpath)

    filenames = ["gene_presence_absence.csv", "Lb_curva_tree.nwk"]
    for filename in filenames:
        filepath = os.path.join(PA.pangenome_dir, filename)
        assert os.path.isfile(filepath) == True
        assert os.stat(filepath).st_size > 0

    remove_prokka_testfiles()

    if os.path.isdir(PA.pangenome_dir):
        shutil.rmtree(PA.pangenome_dir)

    G.create_directory(DC.gff_dir)
    G.create_directory(DC.ffn_dir)
    G.create_directory(DC.fna_dir)

def test_end(config):
    def remove_test_files(config):
        test = config.path
        shutil.rmtree(test)
        tmp_path = os.path.join("/", "home", "pipeline", "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        os.chdir(BASE_PATH)

    remove_test_files(config)

if __name__ == "__main__":
    print(msg)
