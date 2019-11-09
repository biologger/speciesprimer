#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import tarfile
import hashlib
import pytest
import shutil
import json

import urllib.error
import urllib.request

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
from basicfunctions import GeneralFunctions as G
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "blastdb", "tmp")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")
reference_dict = os.path.join(ref_data, "reference_config.json")

def create_mock_archives():
    G.create_directory(tmpdir)
    for i in range(0, 7):
        num = "{:02d}".format(i)
        filestart = "ref_prok_rep_genomes." + num
        mockends = [
            ".nhr", ".nni",
            ".nsi", ".nin",
            ".nog", ".nsq",
            ".nnd", ".nsd"]
        archive = "ref_prok_rep_genomes." + num + ".tar.gz"
        archive_path = os.path.join(tmpdir, archive)
        with tarfile.open(archive_path, "w:gz") as tar:
            for ending in mockends:
                filepath = os.path.join(tmpdir, filestart + ending)
                with open(filepath, "w") as f:
                    f.write("")
                tar.add(filepath)

                if num == "04":
                    os.remove(filepath)

        with open(archive_path + ".md5", "w") as f:
            archmd5 = md5Checksum(archive_path)
            f.write(archmd5 + "  " + archive)

    md5_dir = os.path.join(tmpdir, "md5_files")
    G.create_directory(md5_dir)
    for files in os.listdir(tmpdir):
        if files.endswith("tar.gz.md5"):
            fromfile = os.path.join(tmpdir, files)
            tofilepath = os.path.join(md5_dir, files)
            shutil.copy(fromfile, tofilepath)

    chmdfile = os.path.join(md5_dir, "ref_prok_rep_genomes.05.tar.gz.md5")
    with open(chmdfile, "r") as f:
        for line in f:
            info = line.strip().split("  ")[0][0:-1] + "k  " + line.strip().split("  ")[1]
    with open(chmdfile, "w") as f:
        f.write(info)
    nalfile = os.path.join(tmpdir, "ref_prok_rep_genomes.nal")
    with open(nalfile, "w") as f:
        f.write("")

def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
    return m.hexdigest()

def mocked_urlopen_fail(monkeypatch):
    import urllib.error
    import urllib.request

    def mocked(url, data):
        raise urllib.error.HTTPError(url, 400, "Bad request", None, None)

    monkeypatch.setattr(urllib.request, "urlopen", mocked)

    with pytest.raises(urllib.error.HTTPError):
        pass

def prepare_tmp_db():
    t = os.path.join(testfiles_dir, "tmp_config.json")
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)
    shutil.copy(t, tmp_path)

def remove_tmp_db():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)

def change_tmp_db():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    tmp_dict = {
        "new_run":{
            'modus': "continue", "targets": None, "path": "/primerdesign/test"},
        "email":"biologger@protonmail.com",
        "BLAST_DB":{
            "delete": True, "db": "ref_prok_rep_genomes", 
            "path": "/blastdb/tmp", "test": True}}

    with open(tmp_path, "w") as f:
        f.write(json.dumps(tmp_dict))

def test_getblastdb(monkeypatch):
    create_mock_archives()
    import getblastdb
    prepare_tmp_db()
    change_tmp_db()
    def mocked(url, data):
        raise urllib.error.HTTPError(url, 400, "Bad request", None, None)

    monkeypatch.setattr(urllib.request, "urlopen", mocked)
    getblastdb.get_DB(mode="auto")
    remove_tmp_db()


