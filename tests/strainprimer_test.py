#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import shutil
import pytest
import json
import time
import csv
import multiprocessing
from Bio import SeqIO
from Bio import Entrez
from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G
from basicfunctions import ParallelFunctions as P


msg = (
    """
    Works only in the Docker container!!!
    - Start the container
        sudo docker start {Containername}
    - Start an interactive terminal in the container
        sudo docker exec -it {Containername} bash
    - Start the tests in the container terminal
        cd /
        pytest -vv --cov=/pipeline /tests
    """)


# /tests
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")
testdir = os.path.join("/", "primerdesign", "test")

confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "primerdesign", "test"),
    "probe": False, "exception": [], "minsize": 70, "skip_download": True,
    "customdb": dbpath, "assemblylevel": ["all"], "qc_gene": ["rRNA"],
    "virus": False, "genbank": False, "intermediate": True,
    "nontargetlist": ["Lactobacillus sakei"],
    "evalue": 10, "nuc_identity": 0, "runmode": ["strain"], "strains": []}


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
            args.blastseqs, args.probe, args.virus, args.genbank,
            args.evalue, args.nuc_identity, args.runmode, args.strains)

    config.save_config()

    return config


def test_start():
    if os.path.isdir(testdir):
        shutil.rmtree(testdir)

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
            mockseqs = os.path.join(testfiles_dir, "mocktemplate.seqs")
            records = list(SeqIO.parse(mockseqs, "fasta"))
            SeqIO.write(records, f, "fasta")

        return dbpath

    def create_customblastdb():
        cmd = [
            "makeblastdb", "-in", dbpath, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    if os.path.isdir(os.path.dirname(dbpath)):
        shutil.rmtree(os.path.dirname(dbpath))
    test =  os.path.join("/", "primerdesign", "test")
    if os.path.isdir(test):
        shutil.rmtree(test)

    G.create_directory(os.path.dirname(dbpath))
    dbinputfiles()
    create_customblastdb()


def test_skip_pangenome_analysis(config):
    from speciesprimer import PangenomeAnalysis
    PA = PangenomeAnalysis(config)
    G.create_directory(PA.pangenome_dir)
    fromfile = os.path.join(testfiles_dir, "gene_presence_absence.csv")
    tofile = os.path.join(PA.pangenome_dir, "gene_presence_absence.csv")
    if os.path.isfile(tofile):
        os.remove(tofile)
    shutil.copy(fromfile, tofile)
    exitstat = PA.run_pangenome_analysis()
    assert exitstat == 2


def test_Singleton(config):
    from strainprimer import Singletons
    SI = Singletons(config)

    def prepare_tests():
        if os.path.isdir(SI.ffn_dir):
            shutil.rmtree(SI.ffn_dir)
        new_ffn_dir = os.path.join(testfiles_dir, "ffn_files")
        shutil.copytree(new_ffn_dir, SI.ffn_dir)

    def test_get_singlecopy_genes():
        singleton_count = SI.get_singleton_genes()
        assert singleton_count == 2

    def test_write_fasta():
        locustags = SI.get_sequences_from_ffn()
        SI.get_fasta(locustags)
        filename = os.path.join(SI.blast_dir, "singleton_sequences.fas")
        assert os.path.isfile(filename)

    def test_coregene_extract():
        singletonresdir = os.path.join(
                SI.single_dir, "GCF_003254785v1_20190923")
        if os.path.isdir(singletonresdir):
            shutil.rmtree(singletonresdir)
        SI.singleton_seqs = []
        SI.coregene_extract()
        filename = os.path.join(SI.blast_dir, "singleton_sequences.fas")

        assert os.path.isfile(filename)

    def test_run_singleseqs():
        single_dict = SI.run_singleseqs()
        ref = [
            'GCF_003254785v1_btuD_5', 'GCF_003254785v1_group_3360']
        res = list(single_dict.keys())
        res.sort()
        ref.sort()
        assert res == ref
        return single_dict

    def test_BlastParser(config):
        from strainprimer import SingletonBlastParser
        str_unique = SingletonBlastParser(
            config).run_blastparser(single_dict)
        print("str_unique", str_unique)

    prepare_tests()
    test_get_singlecopy_genes()
    test_write_fasta()
    test_coregene_extract()
    single_dict = test_run_singleseqs()
#    test_BlastParser(config)


def test_end():
    os.chdir(pipe_dir)
    if os.path.isdir(testdir):
        shutil.rmtree(testdir)
