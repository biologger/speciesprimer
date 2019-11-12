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
import pytest

# /tests
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")

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

# functions run in parallel do not count for line coverage
# test these functions here

# primerQC
#MFEprimer_template, MFEprimer_nontarget, MFEprimer_assembly, get_seq_fromDB
# make_templateDB, make_assemblyDB

def initiate_pqc():
    from speciesprimer import PrimerQualityControl
    reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
    with open(reffile) as f:
        for line in f:
            primer3dict = json.loads(line)
    pqc = PrimerQualityControl(config, primer3dict)

if __name__ == "__main__":
    print(msg)
