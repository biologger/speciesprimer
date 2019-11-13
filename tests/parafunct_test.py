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

@pytest.fixture
def blapar(config):
    from speciesprimer import BlastParser
    blapar = BlastParser(config, results="primer")
    return blapar

def prepare_QC_testfiles(pqc, config):
    targetdir = os.path.join(config.path, config.target)
    fna_dir = os.path.join(targetdir, "fna_files")
    qc_file = os.path.join(testfiles_dir, "Lb_curva_qc_sequences_ignoreqc.csv")
    qc_out = os.path.join(pqc.summ_dir, "Lb_curva_qc_sequences.csv")
    G.create_directory(fna_dir)
    G.create_directory(pqc.primer_qc_dir)
    if os.path.isfile(qc_out):
        os.remove(qc_out)
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
        # repeat do not write again
        pqc.make_assemblyDB(primer_qc_list)
        pqc.make_templateDB(primer_qc_list)
        assert os.path.isfile(genDB) == True
        assert os.path.isfile(tempDB) == True
        # remove empty files
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

def primerBLAST(pqc, config):
    from speciesprimer import BlastPrep
    prepare_QC_testfiles(pqc, config)
    G.create_directory(pqc.primerblast_dir)
    prep = BlastPrep(
    pqc.primerblast_dir, pqc.primerlist,
    "primer", pqc.config.blastseqs)
    use_cores, inputseqs = prep.run_blastprep()
    reffile = os.path.join(testfiles_dir, "primer_nontargethits.json")
    tofile = os.path.join(pqc.primerblast_dir, "nontargethits.json")
    shutil.copy(reffile, tofile)

    pqc.call_blastparser.run_blastparser("primer")
    return inputseqs

def prepare_blastdb(config):
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    G.create_directory(tmpdir)
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        dbfile = os.path.join(tmpdir, "primer_customdb.fas")
        with open(dbfile, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbfile

    def create_customblastdb(config, infile):
        cmd = [
            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    dbfile = dbinputfiles()
    create_customblastdb(config, dbfile)

def test_get_seq_from_DB(pqc, config, blapar):
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False
    fasta_seqs = []
    primerBLAST(pqc, config)
    maxsize = 25000
    part = 0
    end = (part+1)*maxsize
    db = config.customdb
    prepare_blastdb(config)
    primerBLAST(pqc, config)
    with open(os.path.join(pqc.primerblast_dir, "nontargethits.json")) as f:
        for line in f:
            nonred_dict = json.loads(line)
    nonreddata = blapar.sort_nontarget_sequences(nonred_dict)
    data = nonreddata[part*maxsize:end]
    data.sort()
    for extract in data:
        fasta = blapar.get_seq_fromDB(extract, db)
        fasta_seqs.append(fasta)

    assert fasta_seqs[0][0] == ">NZ_CP020459.1:1005699-1009917 Lactobacillus sakei strain FAM18311 chromosome, complete genome"
    assert len(fasta_seqs) == 97

def test_end(config):
    def remove_test_files(config):
        test = config.path
        shutil.rmtree(test)
        tmp_path = os.path.join("/", "pipeline", "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.chdir(BASE_PATH)
        assert os.path.isdir(test) == False

    remove_test_files(config)

if __name__ == "__main__":
    print(msg)
