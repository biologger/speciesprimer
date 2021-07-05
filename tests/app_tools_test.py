#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# quick test for tools and dependencies (for Docker image generation)

import os
import time
import shutil
import subprocess
from basicfunctions import GeneralFunctions as G

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
tmpdir = os.path.join("/", "primerdesign", "tmp")

deps = []

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
    """
)

def read_shelloutput(cmd):
    outputlist = []
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, shell=True)

    def check_output():
        while True:
            output = process.stdout.readline().decode().strip()
            if output:
                outputlist.append(output)
            else:
                break
    while process.poll() is None:
        check_output()

    return outputlist

def get_version_infos(cmd):
    toolinfo = read_shelloutput(cmd)
    return toolinfo

def test_versions():
    def blast():
        cmd = " ".join(["blastn", "-version"])
        blastn = get_version_infos(cmd)
        assert blastn[0].split(":")[0] == "blastn"
        print("\nFound blast:")
        print(blastn[0])
        deps.append(blastn[0])

    def fasttree():
        cmd = " ".join(["FastTreeMP"])
        fasttree = get_version_infos(cmd)
        assert ("FastTree Version" in str(fasttree[0])) is True
        version = fasttree[0].split("FastTree Version")[1].split(",")[0]
        print("FastTree" + version)
        deps.append("FastTree" + version)

    def roary():
        cmd1 = "roary -h"
        cmd2 = "roary -version"
        roary = get_version_infos(cmd1)
        version = get_version_infos(cmd2)
        assert (
            roary[0] ==
            'Please cite Roary if you use any of the results it produces:')
        print("Roary", version[0])
        deps.append("Roary " + version[0])

    def prokka():
        prokka = get_version_infos("prokka -version")[0]
        assert prokka.startswith("prokka") is True
        print(prokka)
        deps.append(prokka)

    def prank():
        prank = get_version_infos("prank")[0]
        assert prank.startswith("prank") is True
        version = prank.split(". ")[0]
        print(version)
        deps.append(version)

    def primer3():
        p3 = get_version_infos("primer3_core -h")
        assert p3[0].startswith("Copyright") is True
        for item in p3:
            if "This is primer3" in item:
                version = item.split("This is ")[1]
        print(version)
        deps.append(version)

    def emboss():
        version = get_version_infos("consambig -version")[0]
        assert version.startswith("EMBOSS") is True
        print(version)
        deps.append(version)

    def mfold():
        version = get_version_infos("mfold -h")[0]
        assert version.startswith("mfold") is True
        print(version)
        deps.append(version)

    def mfeprimer():
        version = get_version_infos("MFEprimer.py -v")[0]
        assert version.startswith("MFEprimer") is True
        print(version)
        deps.append(version)

    blast()
    fasttree()
    roary()
    prokka()
    prank()
    primer3()
    emboss()
    mfold()
    mfeprimer()

    with open(os.path.join(BASE_PATH, "tested_versions.md"), "w") as f:
        f.write("### Tools and Versions used for testing\n")
        f.write(time.strftime("%Y/%m/%d", time.localtime()) + " (YY/MM/DD)\n")
        for tool in deps:
            f.write(tool + "\n")

def test_installation():

    def MPprimer_dimer_check():
        os.chdir(os.path.join(todir, "MPprimer"))
        input_file = "Lb_curva_asnS_2_P0"
        output_file = "Lb_curva_asnS_2_P0_dimer_out"
        ref_file = "ref_Lb_curva_asnS_2_P0_dimer_out"
        dimer_cmd = [
            "MPprimer_dimer_check.pl", "-f", input_file, "-d", "3",
            ">", output_file]
        output = read_shelloutput(" ".join(dimer_cmd))

        with open(output_file, "r") as f:
            out = f.readlines()[0].split("\t")
        with open(ref_file, "r") as f:
            ref = f.readlines()[0].split("\t")
        assert out.sort() == ref.sort()
        print(output)

    def mfold():
        os.chdir(os.path.join(todir, "mfold"))
        input_file = "g600_3_P0_PCR"
        output_file = "g600_3_P0_PCR.det"
        ref_file = "ref_g600_3_P0_PCR.det"

        mfold_cmd = [
            "mfold", "SEQ=" + input_file, "NA=DNA", "T=60", "MG_CONC=0.003"]
        output = read_shelloutput(" ".join(mfold_cmd))

        with open(output_file, "r", errors="ignore") as f:
            out = f.readlines()
        with open(ref_file, "r", errors="ignore") as f:
            ref = f.readlines()
        assert out.sort() == ref.sort()
        print(output)

    def MFEprimer():
        print("MFEprimer test")
        os.chdir(os.path.join(todir, "MFEprimer"))
        inputfilepath = "test.rna"
        cmd = ["IndexDB.py", inputfilepath, "-k", "9"]
        output = read_shelloutput(" ".join(cmd))
        print(output)

        dbfiles = [
            ".2bit",
            ".sqlite3.db",
            ".uni"]

        for db in dbfiles:
            assert os.path.isfile(inputfilepath + db) is True

        inputfilepath = "test.rna"
        cmd = ["MFEprimer.py", "-i", "p.fa", "-d", "test.rna"]
        output = read_shelloutput(" ".join(cmd))
        assert output[7] == "8 primer sequences"
        assert output[8] == "Database = test.rna"
        print(output)

    G.create_directory(tmpdir)
    fromdir = os.path.join(testfiles_dir, "tools")
    todir = os.path.join(tmpdir, "tools")
    if not os.path.isdir(todir):
        shutil.copytree(fromdir, todir)

    MPprimer_dimer_check()
    mfold()
    MFEprimer()


def test_end():
    def remove_test_files():
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.chdir(BASE_PATH)
        assert os.path.isdir(tmpdir) is False

    remove_test_files()


if __name__ == "__main__":
    print(msg)