#!/usr/bin/env python3
# -*- coding: utf-8 -*-

msg = (
"""
Works only in the Docker container!!!
- Start the container
    sudo docker start {Containername}
- Start an interactive terminal in the container
    sudo docker exec -it {Containername} bash
- Start the tests in the container terminal
    cd /
    pytest -vv --cov=pipeline /tests/
"""
)

import os
import sys
import json
import shutil
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# /tests
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G
from speciesprimer import PrimerQualityControl

testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")

confargs = {
    "ignore_qc": True, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "primerdesign", "test"),
    "probe": False, "exception": None, "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["tuf"],
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

@pytest.fixture
def pqc(config):
    reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
    with open(reffile) as f:
        for line in f:
            primer3dict = json.loads(line)
    pqc = PrimerQualityControl(config, primer3dict)
    return pqc


# functions run in parallel do not count for line coverage
# test these functions here

# primerQC
#MFEprimer_template, MFEprimer_nontarget, MFEprimer_assembly, get_seq_fromDB
# make_templateDB, make_assemblyDB



def prepare_QC_testfiles(pqc, config):
    targetdir = os.path.join(config.path, config.target)
    fna_dir = os.path.join(targetdir, "fna_files")
    G.create_directory(fna_dir)
    G.create_directory(pqc.primer_qc_dir)
    qc_file = os.path.join(testfiles_dir, "Lb_curva_qc_sequences_ignoreqc.csv")
    qc_out = os.path.join(pqc.summ_dir, "Lb_curva_qc_sequences.csv")
    if os.path.isdir(pqc.summ_dir):
        shutil.rmtree(pqc.summ_dir)
    G.create_directory(pqc.summ_dir)
    shutil.copy(qc_file, qc_out)
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
        newfilename = filename.split(".ffn")[0] + ".fna"
        outpath = os.path.join(fna_dir, newfilename)
        mockfna = SeqRecord(
            Seq("".join(sequences)),
            id="MOCK_" + filename.split("_")[1],
            name="MOCK_" + filename.split("_")[1],
            description="Lactobacillus curvatus")
        with open(outpath, "w") as o:
            SeqIO.write(mockfna, o, "fasta")


def test_make_DBs(pqc, config):
    genDB = os.path.join(pqc.primer_qc_dir, "Lb_curva.genomic")
    tempDB = os.path.join(pqc.primer_qc_dir, "template.sequences") 
    prepare_QC_testfiles(pqc, config)
    if pqc.collect_primer() == 0:
        primer_qc_list = pqc.get_primerinfo(pqc.primerlist, "mfeprimer")
        pqc.make_assemblyDB(primer_qc_list)
        pqc.make_templateDB(primer_qc_list)
        assert os.path.isfile(genDB) == True
        assert os.path.isfile(tempDB) == True        
#        # repeat do not write again       
        pqc.make_assemblyDB(primer_qc_list)
        pqc.make_templateDB(primer_qc_list)
        assert os.path.isfile(genDB) == True
        assert os.path.isfile(tempDB) == True
#        # remove empty files
        os.remove(genDB)
        os.remove(tempDB)
        assert os.path.isfile(genDB) == False
        assert os.path.isfile(tempDB) == False
        dbs = [genDB, tempDB]
        dbfiles = [".2bit", ".sqlite3.db", ".uni"]
        for db in dbs:
            for end in dbfiles:
                if os.path.isfile(db + end):
                    os.remove(db + end)
        qc_out = os.path.join(pqc.summ_dir, "Lb_curva_qc_sequences.csv")
        os.remove(qc_out)
        primer_qc_list = []
        pqc.make_assemblyDB(primer_qc_list)
        pqc.make_templateDB(primer_qc_list)
        assert os.path.isfile(genDB) == False
        assert os.path.isfile(tempDB) == False
        
def test_qc_nottrue(pqc, config):
    confargs["ignore_qc"] == False  
    genDB = os.path.join(pqc.primer_qc_dir, "Lb_curva.genomic")
    tempDB = os.path.join(pqc.primer_qc_dir, "template.sequences") 
    primer_qc_list = []
    pqc.make_assemblyDB(primer_qc_list)
    pqc.make_templateDB(primer_qc_list)
    assert os.path.isfile(genDB) == False
    assert os.path.isfile(tempDB) == False

    os.chdir(pipe_dir)
    shutil.rmtree(config.path)
    if os.path.isdir(pqc.summ_dir):
        shutil.rmtree(pqc.summ_dir)


    

if __name__ == "__main__":
    print(msg)
