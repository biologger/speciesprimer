#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import pytest
import json
from basicfunctions import GeneralFunctions as G
from speciesprimer import Config
from speciesprimer import CLIconf

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
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")

testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")
reference_dict = os.path.join(ref_data, "reference_config.json")

# prompts
start = (
    "Create new config files or start pipeline with previously "
    "generated files?\ntype (n)ew or (s)tart:\n")
species = (
    "Please specify target species (comma separated) "
    "or type help for format examples: \n")
path = (
    "Please specify a path to use as the working directory "
    "or hit return to use the current working directory:\n")
skip_tree = (
    "Skip the core gene alignment and tree for visualization and "
    "troubleshooting.\n(y)es/(n)o, default=(n)\n> ")
offline = (
    "Work offline with local genome assemblies?"
    "\n(y)es/(n)o, default=(n)\n> ")
skip_download = (
    "Skip the download of Genomes from NCBI?"
    "\n(y)es/(n)o, default=(n)\n> ")
assemblylevel = (
    "Limit downloads of Genomes to assembly status\n"
    "options: complete, chromosome, scaffold, contig, "
    "all (comma separated), default=all\n> ")
customdb = (
    "Do you want to use a custom database for blastn?\n"
    "Specifiy the absolute filepath of the custom database "
    "(e.g. '/home/blastdb/nontarget.fasta') or hit return"
    " to skip, default=None\n> ")
blastseqs = (
    "Set the number of sequences per BLAST search. "
    "Decrease the number of sequences if BLAST slows down "
    "due to low memory."
    "\noptions = [100, 500, 1000, 2000, 5000]"
    ", default=1000\n> ")
qc_gene = (
    "Gene(s) (comma separated) for BLAST in the initial quality "
    "control step.\noptions: rRNA, tuf, recA, dnaK, pheS"
    ", default=rRNA\n> ")
exception = (
    "Primer binding to this non-target species is tolerated.\n"
    "Provide a species name or hit return to skip:\n> ")
minsize = ("Minimal Amplicon size\ndefault=70\n> ")
maxsize = ("Maximal amplicon size\ndefault=200\n> ")
probe = (
    "Do you want primer3 to design an internal probe? [Experimental!]"
    "\n(y)es/(n)o, default=(n)\n> ")
mfold = (
    "Delta G threshold for secondary structures in PCR products"
    " at 60 degree Celsius calculated by mfold\ndefault=-3.0\n> ")
mpprimer = (
    "Mpprimer threshold for delta G values calculated by "
    "MPprimer_dimer_check.pl\ndefault=-3.5\n> ")
mfethreshold = (
    "MFEprimer threshold for nontarget sequence PPC "
    "higher values mean more stringent selection.\ndefault=90\n> ")
ignore_qc = (
    "Do you want to include genomes that did not"
    " pass quality control?\ndefault=(n)\n> ")
blastdbv5 = (
    "Do you have the Version 5 of the BLAST DB? \ndefault=(n)\n> ")
intermediate = (
    "Do you want to keep intermediate files?\ndefault=(n)\n> ")
nolist = (
    "Do you want to perform specificity check without the "
    "(non-target) species list (for all sequences in the DB)?"
    "\nNot recommended for nt DB! May be used with a custom DB"
    "\ndefault=(n)\n> ")
forall = (
    "Use this value for all targets?\n(y)es/(n)o, default=(y)\n> ")
targets = (
    "Search for config files for (a)ll or (s)elect targets:\n")


def alldef_input(prompt):
    prompt_dict = {
        start: "n", species: "Lactobacillus curvatus", path: "", skip_tree: "",
        offline: "", skip_download: "", assemblylevel: "", customdb: "",
        blastseqs: "", qc_gene: "", exception: [], minsize: "", maxsize: "",
        probe: "", mfold: "", mpprimer: "", mfethreshold: "", ignore_qc: "",
        blastdbv5: "", intermediate: "", nolist: "", forall: ""}
    val = prompt_dict[prompt]
    return val


def offline_input(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus helveticus",
        path: "/primerdesign/test", skip_tree: "",
        offline: "yes", skip_download: "", assemblylevel: "contig",
        customdb: "",
        blastseqs: "", qc_gene: "", exception: [], minsize: "", maxsize: "",
        probe: "", mfold: "", mpprimer: "", mfethreshold: "", ignore_qc: "",
        blastdbv5: "", intermediate: "", nolist: "", forall: "y"}
    val = prompt_dict[prompt]
    return val


def offline_input2(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus helveticus",
        path: "/primerdesign/test", skip_tree: "", offline: "yes",
        skip_download: "",
        assemblylevel: "contig", customdb: "", blastseqs: "", qc_gene: "",
        exception: [], minsize: "", maxsize: "", probe: "", mfold: "",
        mpprimer: "", mfethreshold: "", ignore_qc: "", blastdbv5: "",
        intermediate: "", nolist: "", forall: "n"}
    val = prompt_dict[prompt]
    return val


def nodef_input(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus helveticus",
        path: "primerdesign/test", skip_tree: "YES", offline: "",
        skip_download: "y", assemblylevel: "all",
        customdb: "/primerdesign/tmp/customdb.fas",
        blastseqs: "2000", qc_gene: "pheS, dnaK",
        exception: "Lactobacillus sunkii",
        minsize: "60", maxsize: "300", probe: "y", mfold: "-2.5",
        mpprimer: "-3.0", mfethreshold: "100", ignore_qc: "y", blastdbv5: "y",
        intermediate: "y", nolist: "y", forall: "n"}
    val = prompt_dict[prompt]
    return val


def wrong_input(prompt):
    prompt_dict = {
        assemblylevel: "alle", customdb: "customdb.fas",
        blastseqs: "66", qc_gene: "tufrrna",
        exception: "Lactobacillus sunkii",
        minsize: "50.36", maxsize: "twohundred", mfold: "minusthree",
        mpprimer: "minusthreepointfive", mfethreshold: "99.9"}
    val = prompt_dict[prompt]
    return val


def wrong_input2(prompt):
    prompt_dict = {
        species: "help", assemblylevel: "alle, complete, contig",
        customdb: "customdb.fas", blastseqs: "forty", qc_gene: "tufrrna",
        exception: "Lactobacillus_sunkii",
        minsize: "50.36", maxsize: "twohundred", mfold: "minusthree",
        mpprimer: "minusthreepointfive", mfethreshold: "99.9"}
    val = prompt_dict[prompt]
    return val


def nodef_input2(prompt):
    prompt_dict = {
        start: "n",
        species: "Lactobacillus curvatus, Lactobacillus helveticus",
        path: "primerdesign/test", skip_tree: "YES", offline: "",
        assemblylevel: "all", customdb: "/primerdesign/tmp/customdb.fas",
        blastseqs: "2000", qc_gene: "pheS, dnaK",
        exception: "Lactobacillus sunkii, Lactobacillus helveticus", minsize: "60", maxsize: "300",
        probe: "y", mfold: "-2.5", mpprimer: "-3.0",
        mfethreshold: "100", ignore_qc: "y", blastdbv5: "y", intermediate: "y",
        skip_download: "y", nolist: "y", forall: "y"}
    val = prompt_dict[prompt]
    return val


def start_input(prompt):
    prompt_dict = {start: "s", targets: "a", path: "/primerdesign/test"}
    val = prompt_dict[prompt]
    return val


def fail_startinput(prompt):
    prompt_dict = {start: "s", targets: "q", path: "/primerdesign/test"}
    val = prompt_dict[prompt]
    return val


def start_wronginput(prompt):
    prompt_dict = {species: "Lactobacillus bifermentans"}
    val = prompt_dict[prompt]
    return val


def save_wronginput(prompt):
    prompt_dict = {
        start: "s", targets: "s", species: "Lactobacillus curvatus",
        path: "/primerdesign/test"}
    val = prompt_dict[prompt]
    return val


def start_oneinput(prompt):
    prompt_dict = {
        start: "s", targets: "s", species: "Lactobacillus curvatus",
        path: "/primerdesign/test"}
    val = prompt_dict[prompt]
    return val


def start_unknown(prompt):
    prompt_dict = {
        start: "s", targets: "s", species: "Lactobacillus sunkii",
        path: "/primerdesign/test"}
    val = prompt_dict[prompt]
    return val


def fail_input(prompt):
    if prompt:
        return "q"


class batchassist_mock():
    def __init__(self, good_input, bad_input):
        self.good_input = good_input
        self.bad_input = bad_input
        self.mocked = []

    def prompt_input(self, prompt):
        if prompt in self.mocked:
            result = self.good_input(prompt)
        else:
            try:
                result = self.bad_input(prompt)
                self.mocked.append(prompt)
            except KeyError:
                result = self.good_input(prompt)
        return result


def compare_configfiles(reference, conffile, key):
    with open(reference) as f:
        for line in f:
            refdict = json.loads(line)
    with open(conffile) as f:
        for line in f:
            confdict = json.loads(line)
    assert refdict[key] == confdict


def get_config_from_file(conf_from_file):
    configfilepaths = []
    targets = conf_from_file.get_targets()
    nontargetlist = []
    targets.sort()
    for target in targets:
        (
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist,
            offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5
                ) = conf_from_file.get_config(target)

        config = CLIconf(
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, nontargetlist, skip_tree,
            nolist, offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5)

        config.save_config()

        cfilepath = os.path.join(
                config.path, config.target, "config", "config.json")
        configfilepaths.append(cfilepath)

    print(configfilepaths)
    return configfilepaths


def test_sys_exit(monkeypatch):
    monkeypatch.setattr('builtins.input', fail_input)
    with pytest.raises(SystemExit):
        conf_from_file = Config()

    monkeypatch.setattr('builtins.input', fail_startinput)
    with pytest.raises(SystemExit):
        conf_from_file = Config()

# mockprimer_csv = os.path.join(path, "Summary", target, abbr + "_primer.csv")


def test_default_input(monkeypatch):
    monkeypatch.setattr('builtins.input', alldef_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "0")
    defdir = os.path.join("/", "Lactobacillus_curvatus")
    if os.path.isdir(defdir):
        shutil.rmtree(defdir)

def test_offline(monkeypatch):
    monkeypatch.setattr('builtins.input', offline_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "1")
    if os.path.isfile(configpath):
        os.remove(configpath)
    configpath = get_config_from_file(conf_from_file)[1]
    compare_configfiles(reference_dict, configpath, "2")
    if os.path.isfile(configpath):
        os.remove(configpath)
    monkeypatch.setattr('builtins.input', offline_input2)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "1")
    if os.path.isfile(configpath):
        os.remove(configpath)
    configpath = get_config_from_file(conf_from_file)[1]
    compare_configfiles(reference_dict, configpath, "2")
    if os.path.isfile(configpath):
        os.remove(configpath)


def test_nodefault(monkeypatch):
    G.create_directory(tmpdir)
    dbpath_tmp = os.path.join(tmpdir, "customdb.fas.nsq")
    with open(dbpath_tmp, "w") as f:
        f.write(">mockDB")
    try:
        monkeypatch.setattr('builtins.input', nodef_input)
        conf_from_file = Config()
        configpath = get_config_from_file(conf_from_file)[0]
        compare_configfiles(reference_dict, configpath, "3")
        if os.path.isfile(configpath):
            os.remove(configpath)
        configpath = get_config_from_file(conf_from_file)[1]
        compare_configfiles(reference_dict, configpath, "4")
        if os.path.isfile(configpath):
            os.remove(configpath)
        monkeypatch.setattr('builtins.input', nodef_input2)
        conf_from_file = Config()
        configpath = get_config_from_file(conf_from_file)[0]
        compare_configfiles(reference_dict, configpath, "8")
        configpath = get_config_from_file(conf_from_file)[1]
        compare_configfiles(reference_dict, configpath, "9")

    finally:
        if os.path.isfile(dbpath_tmp):
            os.remove(dbpath_tmp)    


def test_start_all(monkeypatch):
    monkeypatch.setattr('builtins.input', start_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "8")
    if os.path.isfile(configpath):
        os.remove(configpath)
    configpath = get_config_from_file(conf_from_file)[1]
    compare_configfiles(reference_dict, configpath, "9")
    if os.path.isfile(configpath):
        os.remove(configpath)


def test_start_unknown(monkeypatch):
    print("Start unknown")
    mock = batchassist_mock(start_oneinput, start_unknown)
    monkeypatch.setattr('builtins.input', mock.prompt_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "8")


def test_start_wrong(monkeypatch):
    print("Start wrong")
    mock = batchassist_mock(save_wronginput, start_wronginput)
    monkeypatch.setattr('builtins.input', mock.prompt_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "8")
    if os.path.isfile(configpath):
        os.remove(configpath)


def test_not_valid_input(monkeypatch):
    print("test not valid input")
    dbpath_tmp = os.path.join(tmpdir, "customdb.fas.nsq")
    with open(dbpath_tmp, "w") as f:
        f.write(">mockDB")
    try:
        mock = batchassist_mock(nodef_input, wrong_input)
        monkeypatch.setattr('builtins.input', mock.prompt_input)
        conf_from_file = Config()
        configpath = get_config_from_file(conf_from_file)[0]
        compare_configfiles(reference_dict, configpath, "3")
        if os.path.isfile(configpath):
            os.remove(configpath)
        configpath = get_config_from_file(conf_from_file)[1]
        compare_configfiles(reference_dict, configpath, "4")
        if os.path.isfile(configpath):
            os.remove(configpath)
        mock = batchassist_mock(nodef_input, wrong_input2)
        monkeypatch.setattr('builtins.input', mock.prompt_input)
        conf_from_file = Config()
        configpath = get_config_from_file(conf_from_file)[0]
        compare_configfiles(reference_dict, configpath, "7")
        if os.path.isfile(configpath):
            os.remove(configpath)
    finally:
        if os.path.isfile(dbpath_tmp):
            os.remove(dbpath_tmp)


def test_incomplete_config(monkeypatch):
    testfile = os.path.join(
        "/", "primerdesign", "test", "Lactobacillus_curvatus",
        "config", "config.json")
    incomplete_dict = {
            "target": "lactobacillus_curvatus", "path": "/primerdesign/test"}
    with open(testfile, "w") as f:
        f.write(json.dumps(incomplete_dict))
    monkeypatch.setattr('builtins.input', start_input)
    conf_from_file = Config()
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "5")


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
        "email": "biologger@protonmail.com"}
    with open(tmp_path, "w") as f:
        f.write(json.dumps(tmp_dict))


def change_db_again():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    tmp_dict = {
        "new_run": {
            'modus': "continue",
            "targets": ["Lactobacillus_curvatus", "Lactobacillus_sunkii"],
            "path": "/primerdesign/test"},
        "email": "biologger@protonmail.com"}
    with open(tmp_path, "w") as f:
        f.write(json.dumps(tmp_dict))


def test_autorun():
    from speciesprimer import auto_run
    prepare_tmp_db()
    targets, conf_from_file, use_configfile = auto_run()
    assert targets == ["Lactobacillus_curvatus"]
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "6")
    change_tmp_db()
    targets, conf_from_file, use_configfile = auto_run()
    targets.sort()
    assert targets == ["Lactobacillus_curvatus", "Lactobacillus_helveticus"]
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "6")
    change_db_again()
    targets, conf_from_file, use_configfile = auto_run()
    assert targets == ["Lactobacillus_curvatus"]
    configpath = get_config_from_file(conf_from_file)[0]
    compare_configfiles(reference_dict, configpath, "6")
    remove_tmp_db()


def test_shellrun(monkeypatch):
    conf1 = os.path.join(
        "/", "primerdesign", "test", "Lactobacillus_curvatus",
        "config", "config.json")
    conf2 = os.path.join(
        "/", "primerdesign", "test", "Lactobacillus_helveticus",
        "config", "config.json")
    test = os.path.join("/", "primerdesign", "test")
    if os.path.isdir(test):
        shutil.rmtree(test)
    monkeypatch.setattr('builtins.input', offline_input)
    from speciesprimer import main
    main()
    assert os.path.isfile(conf1) is True
    assert os.path.isfile(conf2) is True


def test_end():
    def remove_test_files():
        test = os.path.join("/", "primerdesign", "test")
        shutil.rmtree(test)
        tmp_path = os.path.join("/", "pipeline", "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.chdir(BASE_PATH)
        assert os.path.isdir(test) is False

    remove_test_files()


if __name__ == "__main__":
    print(msg)
