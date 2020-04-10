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
        pytest -vv --cov=pipeline /tests/
    """)


# /tests
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")

confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "primerdesign", "test"),
    "probe": False, "exception": [], "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["rRNA"],
    "virus": False, "genbank": False, "intermediate": True,
    "nontargetlist": ["Lactobacillus sakei"],
    "evalue": 500, "nuc_identity": 0, "runmode": ["species"], "strains": []}


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


def generate_tmp_config(key, value, key2, value2):
    from speciesprimer import CLIconf
    confargs[key] = value
    confargs[key2] = value2
    args = AttrDict(confargs)
    nontargetlist = []
    config = CLIconf(
            args.minsize, args.maxsize, args.mpprimer, args.exception,
            args.target, args.path, args.intermediate,
            args.qc_gene, args.mfold, args.skip_download,
            args.assemblylevel, nontargetlist,
            args.skip_tree, args.nolist, args.offline,
            args.ignore_qc, args.mfethreshold, args.customdb,
            args.blastseqs, args.probe, args.virus, args.genbank,
            args.evalue, args.nuc_identity, args.runmode, args.strains)

    return config

def test_filter_nonreddict(config):
    from speciesprimer import BlastParser
    nonred = os.path.join(testfiles_dir, "nontargethits.json")
    with open(nonred) as f:
        for line in f:
            nonred_dict = json.loads(line)

    ref_diff = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 4, 4, 26, 5, 5, 2, 0, 4, 1, 16, 0, 0, 1, 0, 3, 5, 13, 9, 0, 7,
        5, 7, 3, 4, 0, 3, 0, 7, 0, 7, 57, 10, 0, 4, 1, 0],
        [0, 0, 0, 4, 4, 26, 5, 5, 2, 0, 4, 1, 16, 0, 0, 1, 0, 3, 5, 13, 9, 0, 7,
        5, 7, 3, 4, 0, 3, 49, 7, 0, 7, 57, 10, 0, 4, 1, 4]]

    testsett = [
        [500, 0], [10, 0], [10, 75], [10, 80], [10, 90], [1e-20, 100]]
    print("config evalue, nucleotide identity")
    for i, sett in enumerate(testsett):

        config = generate_tmp_config(
                                "evalue", sett[0], "nuc_identity", sett[1])
        blapa = BlastParser(config)
        print(blapa.config.evalue, blapa.config.nuc_identity)

        filtereddict = blapa.filter_nonreddict(nonred_dict)

        diffs = []
        for k, v in nonred_dict.items():
            diffs.append(len(v) - len(filtereddict[k]))

        assert diffs.sort() == ref_diff[i].sort()


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
        assert os.path.isdir(test) is False

    remove_test_files(config)

if __name__ == "__main__":
    print(msg)
