#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import tarfile
import hashlib
import pytest
import shutil
import json
import urllib.error
import urllib.request
import getblastdb
from getblastdb import config
from basicfunctions import GeneralFunctions as G


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

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "blastdb", "tmp")
downdir = os.path.join(tmpdir, "mockfiles", "download")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")
reference_dict = os.path.join(ref_data, "reference_config.json")

@pytest.fixture
def delete_config():
    db = "ref_prok_rep_genomes"
    delete = True
    test = True
    down_c = config(db, tmpdir, delete, test)
    return down_c

@pytest.fixture
def keep_config():
    db = "ref_prok_rep_genomes"
    delete = False
    test = True
    down_c = config(db, tmpdir, delete, test)
    return down_c

def create_mock_archives():
    G.create_directory(tmpdir)
    G.create_directory(downdir)
    refhtml = os.path.join(testfiles_dir, "download.html")
    htmlfile = os.path.join(tmpdir, "mockfiles","download.html")
    shutil.copy(refhtml, htmlfile)
    with open(os.path.join(tmpdir, "ref_prok_rep_genomes.html.tmp"), "w") as f:
        f.write("mock")
    for i in range(0, 6):
        num = "{:02d}".format(i)
        filestart = "ref_prok_rep_genomes." + num
        mockends = [
            ".nhr", ".nni",
            ".nsi", ".nin",
            ".nog", ".nsq",
            ".nnd", ".nsd"]
        archive = "ref_prok_rep_genomes." + num + ".tar.gz"
        os.chdir(tmpdir)
        with tarfile.open(archive, "w:gz") as tar:
            for ending in mockends:
                filepath = filestart + ending
                with open(filepath, "w") as f:
                    f.write("")
                tar.add(filepath)

        with open(archive + ".md5", "w") as f:
            archmd5 = md5Checksum(archive)
            f.write(archmd5 + "  " + archive)

        archfile = "ref_prok_rep_genomes.03.tar.gz"
        tmpfile = os.path.join(downdir, "ref_prok_rep_genomes.03.tar.gz")
        if os.path.isfile(archfile):
            shutil.move(archfile, tmpfile)

        archfile = "ref_prok_rep_genomes.05.tar.gz"
        tmpfile = os.path.join(downdir, "ref_prok_rep_genomes.05.tar.gz")
        if os.path.isfile(archfile):
            shutil.copy(archfile, tmpfile)

        for files in os.listdir(tmpdir):
            if (
                files.startswith("ref_prok_rep_genomes.03")
                and not files.endswith(".md5")
            ):
                os.remove(files)

    md5_dir = os.path.join(tmpdir, "md5_files")
    G.create_directory(md5_dir)
    for files in os.listdir(tmpdir):
        if files.endswith("tar.gz.md5"):
            fromfile = os.path.join(tmpdir, files)
            tofilepath = os.path.join(md5_dir, files)
            downfiles = os.path.join(downdir, files)
            shutil.copy(fromfile, tofilepath)
            shutil.copy(fromfile, downfiles)

    chmdfile = os.path.join(md5_dir, "ref_prok_rep_genomes.05.tar.gz.md5")
    with open(chmdfile, "r") as f:
        for line in f:
            info = (
                line.strip().split("  ")[0][0:-1] + "k  " +
                line.strip().split("  ")[1])
    with open(chmdfile, "w") as f:
        f.write(info)
    nalfile = os.path.join(tmpdir, "ref_prok_rep_genomes.nal")
    with open(nalfile, "w") as f:
        f.write("")
    os.remove("ref_prok_rep_genomes.04.tar.gz")


def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
    return m.hexdigest()


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
        "new_run": {
            'modus': "continue", "targets": None,
            "path": "/primerdesign/test"},
        "email": "biologger@protonmail.com",
        "BLAST_DB": {
            "delete": True, "db": "ref_prok_rep_genomes",
            "path": "/blastdb/tmp", "test": True}}

    with open(tmp_path, "w") as f:
        f.write(json.dumps(tmp_dict))


def check_all_files():
    for i in range(0, 6):
        num = "{:02d}".format(i)
        filestart = "ref_prok_rep_genomes." + num
        mockends = [
            ".nhr", ".nni",
            ".nsi", ".nin",
            ".nog", ".nsq",
            ".nnd", ".nsd"]
        for ending in mockends:
            filepath = filestart + ending
            assert os.path.isfile(filepath) is True


def check_part_five(outcome=True):
    num = "05"
    filestart = "ref_prok_rep_genomes." + num
    mockends = [
        ".nhr", ".nni",
        ".nsi", ".nin",
        ".nog", ".nsq",
        ".nnd", ".nsd"]
    for ending in mockends:
        filepath = filestart + ending
        assert os.path.isfile(filepath) == outcome


def test_commandline():
    parser = getblastdb.commandline()
    args = parser.parse_args([])
    # default
    assert args.database == 'nt'
    assert args.dbpath is None
    assert args.delete is False
    args = parser.parse_args([
            "--database", "ref_prok_rep_genomes", "--dbpath", "/blastdb",
            "--delete"])
    assert args.database == "ref_prok_rep_genomes"
    assert args.dbpath == "/blastdb"
    assert args.delete is True
    args = parser.parse_args([
            "--database", "ref_prok_rep_genomes", "--dbpath", "/blastdb"])
    assert args.database == "ref_prok_rep_genomes"
    assert args.dbpath == "/blastdb"
    assert args.delete is False
    with pytest.raises(SystemExit):
        args = parser.parse_args([
                "database", "ref_prok_rep_genomes", "dbpath", "/blastdb"])


def test_run():
    try:
        remfile = os.path.join(
                tmpdir, "ref_prok_rep_genomes.05.tar.gz.md5")
        create_mock_archives()
        cmd = [
                "getblastdb.py", "--database", "ref_prok_rep_genomes",
                "--dbpath", "/blastdb/tmp", "--test", "--delete"]
        os.chdir(tmpdir)
        os.remove(remfile)

        for files in os.listdir(tmpdir):
            if (
                files.startswith("ref_prok_rep_genomes.05")
                and not files.endswith(".md5")
            ):
                os.remove(files)

        G.run_subprocess(cmd)

        file1 = os.path.join(
                tmpdir, "md5_files", "ref_prok_rep_genomes.00.tar.gz.md5")
        file2 = os.path.join(
                tmpdir, "md5_files", "ref_prok_rep_genomes.05.tar.gz.md5")
        assert os.path.isfile(file1) is True
        assert os.path.isfile(file2) is True
        check_all_files()
    finally:
        os.chdir("..")
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)

def test_getmd5files(delete_config):
    create_mock_archives()
    getblastdb.get_md5files(delete_config)
    file1 = os.path.join(tmpdir, "ref_prok_rep_genomes.00.tar.gz.md5")
    file2 = os.path.join(tmpdir, "ref_prok_rep_genomes.05.tar.gz.md5")
    assert os.path.isfile(file1) is True
    assert os.path.isfile(file2) is True
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)

def test_skipdownload(delete_config):
    create_mock_archives()
    os.chdir(tmpdir)
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    with open(chmdfile, "r") as f:
        for line in f:
            info = (
                line.strip().split("  ")[0][0:-1] + "k  " +
                line.strip().split("  ")[1])
    with open(chmdfile, "w") as f:
        f.write(info)
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_skipdownload_nomd5(delete_config):
    create_mock_archives()
    os.chdir(tmpdir)
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    refmd5 = os.path.join("md5_files", chmdfile)
    if os.path.isfile(refmd5):
        os.remove(refmd5)
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_archiveextract(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and not files.endswith(".md5")
        ):
            os.remove(files)
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    assert os.path.isfile("ref_prok_rep_genomes.05.tar.gz") is False
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_onlymd5(keep_config):

    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    refmd5 = os.path.join("md5_files", chmdfile)
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and not files.endswith(".md5")
        ):
            os.remove(files)
    shutil.copy(chmdfile, refmd5)
    getblastdb.download_from_ftp([chmdfile], keep_config)
    check_part_five(outcome=True)
    assert os.path.isfile("ref_prok_rep_genomes.05.tar.gz") is True
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_no_old_md5(keep_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and not files.endswith(".md5")
        ):
            os.remove(files)
    shutil.rmtree("md5_files")
    G.create_directory("md5_files")
    getblastdb.download_from_ftp([chmdfile], keep_config)
    check_part_five(outcome=True)
    assert os.path.isfile("ref_prok_rep_genomes.05.tar.gz") is True
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_checksum_incorrect(delete_config):
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and not files.endswith(".md5")
        ):
            os.remove(files)
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    with open(chmdfile, "r") as f:
        for line in f:
            info = (
                line.strip().split("  ")[0][0:-1] + "q67  " +
                line.strip().split("  ")[1])
    with open(chmdfile, "w") as f:
        f.write(info)
    with pytest.raises(Exception):
        getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=False)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_compare_md5_archive(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and "tar.gz" not in files
        ):
            os.remove(files)
    shutil.rmtree("md5_files")
    G.create_directory("md5_files")
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_compare_md5_files_witharchive(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    create_mock_archives()
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and "tar.gz" not in files
        ):
            os.remove(files)
    archive = os.path.join(downdir, "ref_prok_rep_genomes.05.tar.gz")
    shutil.copy(archive, "ref_prok_rep_genomes.05.tar.gz")
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_compare_md5_archive_changed(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    archive = chmdfile.split(".md5")[0]
    create_mock_archives()
    shutil.rmtree("md5_files")
    G.create_directory("md5_files")
    os.chdir(tmpdir)
    with tarfile.open(archive, "w:gz") as tar:
        for files in os.listdir(tmpdir):
            if files.startswith("ref_prok_rep_genomes.05"):
                tar.add(files)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and "tar.gz" not in files
        ):
            os.remove(files)
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_compare_corrupted_md5(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    create_mock_archives()
    shutil.rmtree("md5_files")
    G.create_directory("md5_files")
    os.chdir(tmpdir)
    for files in os.listdir(tmpdir):
        if (
            files.startswith("ref_prok_rep_genomes.05")
            and "tar.gz" not in files
        ):
            os.remove(files)
    with open(chmdfile, "r") as f:
        for line in f:
            info = (
                line.strip().split("  ")[0][0:-1] + "q67  " +
                line.strip().split("  ")[1])
    with open(chmdfile, "w") as f:
        f.write(info)

    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_compare_md5_files(delete_config):
    chmdfile = "ref_prok_rep_genomes.05.tar.gz.md5"
    refmd5 = os.path.join("md5_files", chmdfile)
    create_mock_archives()
    os.chdir(tmpdir)
    with open(refmd5, "r") as f:
        for line in f:
            info = (
                line.strip().split("  ")[0][0:-1] + "q67  " +
                line.strip().split("  ")[1])
    with open(refmd5, "w") as f:
        f.write(info)
    getblastdb.download_from_ftp([chmdfile], delete_config)
    check_part_five(outcome=True)
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_getblastdb_auto():
    create_mock_archives()
    prepare_tmp_db()
    change_tmp_db()
    getblastdb.get_DB(mode="auto")
    check_all_files()
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_badrequest(monkeypatch):
    def mocked(url, data):
        raise urllib.error.HTTPError(url, 400, "Bad request", None, None)
    create_mock_archives()
    monkeypatch.setattr(urllib.request, "urlopen", mocked)
    with pytest.raises(Exception):
        getblastdb.get_DB(mode="auto")
    num = "03"
    filestart = "ref_prok_rep_genomes." + num
    mockends = [
        ".nhr", ".nni",
        ".nsi", ".nin",
        ".nog", ".nsq",
        ".nnd", ".nsd"]
    for ending in mockends:
        filepath = filestart + ending
        assert os.path.isfile(filepath) is False
    remove_tmp_db()
    os.chdir("..")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    assert os.path.isdir(tmpdir) is False


if __name__ == "__main__":
    print(msg)
