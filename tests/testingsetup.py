#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pytest
from basicfunctions import HelperFunctions as H

msg = (
    """
    Works only in the Docker container!!!
    - Start the container
        sudo docker start {Containername}
    - Start an interactive terminal in the container
        sudo docker exec -u 0 -it {Containername} bash
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


class MockingBird():
    def __init__(self):
        self.count = 0

    def mock_syn(self, *args, **kwargs):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock01.xml")
        f = open(mockfile, 'rb')
        return f

    def mock_getsummary(self, *args, **kwargs):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "getsummarymock.xml")
        f = open(mockfile, 'rb')
        return f

    def mock_getlinks(self, *args, **kwargs):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "getlinksmock.xml")
        f = open(mockfile, 'rb')
        return f

    def run(self, *args, **kwargs):
        if self.count == 0:
            self.count += 1
            return self.mock_syn(*args)

        if self.count == 1:
            self.count += 1
            return  self.mock_getsummary(*args, **kwargs)

        if self.count == 2:
            self.count += 1
            return  self.mock_getlinks(*args, **kwargs)


@pytest.fixture
def config():
    from speciesprimer import CLIconf
    args = AttrDict(confargs)
    nontargetlist = H.create_non_target_list(args.target)
    config = CLIconf(
            args.minsize, args.maxsize, args.mpprimer, args.exception,
            args.target, args.path, args.intermediate,
            args.qc_gene, args.mfold, args.skip_download,
            args.assemblylevel, 
            args.skip_tree, args.nolist, args.offline,
            args.ignore_qc, args.mfethreshold, args.customdb,
            args.blastseqs, args.probe, args.virus, args.genbank,
            args.evalue, args.nuc_identity, args.runmode, args.strains,
            nontargetlist)

    config.save_config()

    return config