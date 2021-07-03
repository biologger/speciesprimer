#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import pytest
import json
from Bio import Entrez
from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G
import filecmp

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

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def test_BLASTDB_check():
    confargs = {"blastdbv5": False, "customdb": "/blastdb/not_a_DB.fas"}
    config = AttrDict(confargs)
    with pytest.raises(Exception):
        H.BLASTDB_check(config)


def test_advancedconfig_fromfile():
    CONFFILE = os.path.join(testfiles_dir, "adconfig", "advanced_config.json")
    H.advanced_pipe_config(CONFFILE)
    filenames = [
        "no_blast.gi", "genus_abbrev.csv", "p3parameters", "species_list.txt"]

    for filename in filenames:
        ref = os.path.join(ref_data, "advanced_config", filename)
        newfile = os.path.join(dict_path, filename)
        assert filecmp.cmp(newfile, ref) is True
        deffile =  os.path.join(dict_path, "default", filename)
        os.remove(newfile)
        shutil.copy(deffile, newfile)

    G.create_directory(tmpdir)
    for filename in filenames:
        ref = os.path.join(ref_data, "advanced_config", filename)
        tmp = os.path.join(tmpdir, filename)
        shutil.copy(ref, tmp)
    CONFFILE = os.path.join(testfiles_dir, "adconfig", "advanced_config2.json")
    H.advanced_pipe_config(CONFFILE)
    for filename in filenames:
        ref = os.path.join(ref_data, "advanced_config", filename)
        newfile = os.path.join(dict_path, filename)
        assert filecmp.cmp(newfile, ref) is True
        deffile =  os.path.join(dict_path, "default", filename)
        os.remove(newfile)
        shutil.copy(deffile, newfile)

    CONFFILE = os.path.join(testfiles_dir, "adconfig","advanced_config3.json")
    tmp = os.path.join(tmpdir, "p3parameters")
    wrong_ext = os.path.join(tmpdir, "p3parameters.txt")
    shutil.copy(tmp, wrong_ext)
    exitstat = H.advanced_pipe_config(CONFFILE)
    assert exitstat == 1

    CONFFILE = os.path.join(testfiles_dir, "adconfig", "advanced_config4.json")
    exitstat = H.advanced_pipe_config(CONFFILE)
    assert exitstat == 1
    for filename in filenames:
        ref = os.path.join(ref_data, "advanced_config", filename)
        newfile = os.path.join(dict_path, filename)
        deffile =  os.path.join(dict_path, "default", filename)
        os.remove(newfile)
        shutil.copy(deffile, newfile)

    CERTFILE = os.path.join(testfiles_dir, "adconfig", "mock_certificate.crt")
    CERT_REF = os.path.join(tmpdir, "mock_certificate.crt")
    CONFFILE = os.path.join(tmpdir, "adconfig.json")
    CERT_EXT = os.path.join(tmpdir, "mock_certificate.txt")
    shutil.copy(CERTFILE, CERT_REF)
    shutil.copy(CERT_REF, CERT_EXT)
    conf = {"certificate": CERT_REF}
    with open(CONFFILE, "w") as f:
        f.write(json.dumps(conf))

    exitstat = H.advanced_pipe_config(CONFFILE)
    assert exitstat == 0

    os.remove(CONFFILE)

    conf = {"certificate": CERT_EXT}
    with open(CONFFILE, "w") as f:
        f.write(json.dumps(conf))
    exitstat = H.advanced_pipe_config(CONFFILE)
    assert exitstat == 1

    os.remove(CONFFILE)

    conf = {"certificate": tmpdir}
    with open(CONFFILE, "w") as f:
        f.write(json.dumps(conf))
    exitstat = H.advanced_pipe_config(CONFFILE)
    assert exitstat == 1

    shutil.rmtree(tmpdir)


def prepare_tmp_db():
    t = os.path.join(BASE_PATH, "testfiles", "tmp_config.json")
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)
    shutil.copy(t, tmp_path)


def change_tmp_db():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    tmp_dict = {'modus': "continue", "targets": None}
    with open(tmp_path, "w") as f:
        f.write(json.dumps(tmp_dict))


def remove_tmp_db():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    if os.path.isfile(tmp_path):
        os.remove(tmp_path)


def read_email_from_db():
    tmp_path = os.path.join(pipe_dir, "tmp_config.json")
    if os.path.isfile(tmp_path):
        with open(tmp_path) as f:
            for line in f:
                tmp_db = json.loads(line)
        try:
            mail = tmp_db['email']
            if '@' and '.' in mail:
                email = mail.strip()
                return email
        except KeyError:
            return None
    else:
        return False


class Entrezmail():
    def __init__(self, outcome):
        self.repeat = 0
        self.outcome = outcome

    def email_input(self, prompt):
        user_email = (
                "To make use of NCBI's E-utilities, "
                "Please enter your email address. \n")
        if self.repeat == 0:
            prompt_dict = {user_email: "hfkjdsalgfhlgh.com"}
            assert read_email_from_db() == self.outcome
            self.repeat += 1
        elif self.repeat < 2:
            prompt_dict = {user_email: "hfkjdsa@lgfhlghcom"}
            assert read_email_from_db() == self.outcome
            self.repeat += 1
        else:
            prompt_dict = {user_email: "biologger@protonmail.com"}
            assert read_email_from_db() == self.outcome
        val = prompt_dict[prompt]
        return val


def test_Entrez_email(monkeypatch):
    remove_tmp_db()
    monkeypatch.setattr('builtins.input', Entrezmail(False).email_input)
    H.get_email_for_Entrez(email=None)
    assert read_email_from_db() == "biologger@protonmail.com"
    remove_tmp_db()
    H.get_email_for_Entrez(email="biologgerprton.com")
    assert read_email_from_db() == "biologger@protonmail.com"
    remove_tmp_db()
    H.get_email_for_Entrez(email="biologger@protonmail.com")
    assert read_email_from_db() == "biologger@protonmail.com"
    change_tmp_db()
    monkeypatch.setattr('builtins.input', Entrezmail(None).email_input)
    H.get_email_for_Entrez(email=None)
    assert read_email_from_db() == "biologger@protonmail.com"
    change_tmp_db()
    H.get_email_for_Entrez(email="biologgerprton.com")
    assert read_email_from_db() == "biologger@protonmail.com"
    change_tmp_db()
    H.get_email_for_Entrez(email="biologger@protonmail.com")
    assert read_email_from_db() == "biologger@protonmail.com"
    remove_tmp_db()


def test_subsp_abbrev():
    target = "Lactococcus_lactis_subsp_lactis"
    name = H.abbrev(target)
    assert name == "Lc_lacti_lacti"
    target = "Salmonella_enterica_subsp_enterica"
    name = H.abbrev(target)
    assert name == "Salmo_enter_enter"
    target = "Salmonella_enterica"
    name = H.abbrev(target)
    assert name == "Salmo_enter"


def a_function(inputitem):
    if inputitem == "lactis subsp. lactis":
        raise Exception
    return inputitem


def test_run_parallel_exceptions():
    inputlist = [
            "lactis_subsp_lactis", "lactis subsp. lactis",
            "enterica", "enterica"]
    result = G.run_parallel(a_function, inputlist, args=False, verbosity="")
    ref = ['lactis_subsp_lactis', 'enterica', 'enterica']
    assert result.sort() == ref.sort()


def test_subspecies_handler():
    outcome = [
            "lactis_subsp_lactis", "lactis subsp. lactis",
            "enterica", "enterica",
            "syringae_pv_syringae", "syringae pv. syringae"]
    targets = [
        "Lactococcus_lactis_subsp_lactis", "Salmonella_enterica",
        "Pseudomonas_syringae_pv_syringae"]
    modes = ["underscore", "space"]
    outnum = 0
    for target in targets:
        for mode in modes:
            species = H.subspecies_handler(target, mode=mode)
            assert species == outcome[outnum]
            outnum += 1


def test_check_input(monkeypatch):
    def mock_taxid(db, term):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "esearchmock01.xml")
        f = open(mockfile, 'rb')
        return f
    target = "Lactobacillus_curvatus"
    email = "biologger@protonmail.com"
    monkeypatch.setattr(Entrez, "esearch", mock_taxid)
    taxid = H.check_input(target, email)
    assert taxid == str(28038)


def test_check_input_fail(monkeypatch):
    def mock_taxid(db, term):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "esearchmock02.xml")
        f = open(mockfile, 'rb')
        return f
    target = "Lactobacious_curvatus"
    email = "biologger@protonmail.com"
    monkeypatch.setattr(Entrez, "esearch", mock_taxid)
    taxid = H.check_input(target, email)
    assert taxid is None


def test_check_input_OSError(monkeypatch):
    import urllib.error
    import urllib.request
    def mock_error(db, term):
        raise urllib.error.HTTPError(
                    "Entrez", 400, "Bad request", None, None)
    target = "Lactobacious_curvatus"
    email = "biologger@protonmail.com"
    monkeypatch.setattr(Entrez, "esearch", mock_error)
    with pytest.raises(OSError):
        H.check_input(target, email)


def test_check_species_syn(monkeypatch):
    def mock_syn(db, id):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock01.xml")
        f = open(mockfile, 'rb')
        return f
    target = "Lactobacillus_curvatus"
    email = "biologger@protonmail.com"
    taxid = str(28038)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn == [
        'Bacterium curvatum', 'Lactobacillus curvatus subsp. curvatus',
        'Lactobacillus sp. N55', 'Lactobacillus sp. N61']


def test_check_species_syn_target_sc(monkeypatch):
    def mock_syn(db, id):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock04.xml")
        f = open(mockfile, 'rb')
        return f
    target = "Chlamydophila pneumoniae"
    email = "biologger@protonmail.com"
    taxid = str(28038)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn == ['Chlamydia pneumoniae']
    target = 'Chlamydia pneumoniae'
    email = "biologger@protonmail.com"
    taxid = str(28038)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn == ["Chlamydophila pneumoniae"]


def test_check_species_syn_none(monkeypatch):
    def mock_syn(db, id):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock02.xml")
        f = open(mockfile, 'rb')
        return f
    # not relevant since return is mocked
    target = "Lactobacillus_curvatus"
    email = "biologger@protonmail.com"
    taxid = str(28038)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn is None


def test_check_synonyms_nosyn(monkeypatch):
    def mock_syn(db, id):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock03.xml")
        f = open(mockfile, 'rb')
        return f
    # not relevant since return is mocked
    target = "Lactobacillus_wasatchensis"
    email = "biologger@protonmail.com"
    taxid = str(1335616)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn == ['Lactobacillus sp. WDC04']


def test_check_synonyms_nonesyn(monkeypatch):
    def mock_syn(db, id):
        mockfile = os.path.join(
                testfiles_dir, "entrezmocks", "efetchmock05.xml")
        f = open(mockfile, 'rb')
        return f
    # not relevant since return is mocked
    target = "Lactobacillus_curvatus"
    email = "biologger@protonmail.com"
    taxid = str(28038)
    monkeypatch.setattr(Entrez, "efetch", mock_syn)
    syn = H.check_species_syn(taxid, email, target)
    assert syn is None


def test_create_non_target_list():
    targetdef = []
    target = "Lactococcus_lactis_subsp_lactis"
    nontargetlist = H.create_non_target_list(target)
    for item in nontargetlist:
        if (
            item == "Lactococcus lactis subsp. lactis" or
            item == "Lactococcus lactis subsp lactis"
        ):
            targetdef.append(item)
    assert targetdef == []


def test_rollback():
    testdir = os.path.join("/", "primerdesign", "test")
    target_dir = os.path.join(testdir, "Lactobacillus_curvatus")
    pangenome_dir = os.path.join(target_dir, "Pangenome")
    results_dir = os.path.join(pangenome_dir, "results")
    alignments_dir = os.path.join(results_dir, "alignments")
    consensus_dir = os.path.join(results_dir, "consensus")
    primer_dir = os.path.join(results_dir, "primer")
    primer_qc_dir = os.path.join(primer_dir, "primer_QC")
    blast_dir = os.path.join(results_dir, "blast")
    primerblast_dir = os.path.join(primer_dir, "primerblast")
    G.create_directory(primerblast_dir)
    G.create_directory(alignments_dir)
    G.create_directory(consensus_dir)
    G.create_directory(blast_dir)
    G.create_directory(primer_qc_dir)

    # annotation
    annotdir = "GCF_902362325v1_20191114"
    G.create_directory(os.path.join(target_dir, annotdir))
    outdir = annotdir
    assert os.path.isdir(os.path.join(target_dir, annotdir)) is True
    G.keyexit_rollback("annotation", dp=os.path.join(target_dir, outdir))
    assert os.path.isdir(os.path.join(target_dir, annotdir)) is False
    # p3
    p3_file = os.path.join(primer_dir, "primer3_output")
    output_file = p3_file
    with open(p3_file, "w") as f:
        f.write("PRIMER3MOCK")
    assert os.path.isfile(p3_file) is True
    G.keyexit_rollback("primer3 run", fp=output_file)
    assert os.path.isfile(p3_file) is False
    # alignments
    run_file = os.path.join(results_dir, "run_prank")
    with open(run_file, "w") as f:
        f.write("PRANK MSA MOCK")
    assert os.path.isfile(run_file) is True
    G.keyexit_rollback("Prank MSA run", dp=alignments_dir, fp=run_file)
    assert os.path.isfile(run_file) is False
    assert os.path.isdir(alignments_dir) is False
    # consensus
    run_file = os.path.join(results_dir, "run_consensus")
    with open(run_file, "w") as f:
        f.write("CONSAMBIG MOCK")
    assert os.path.isfile(run_file) is True
    G.keyexit_rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) is False
    assert os.path.isdir(consensus_dir) is False
    # dp only
    run_file = os.path.join(results_dir, "run_consensus")
    G.create_directory(consensus_dir)
    assert os.path.isfile(run_file) is False
    G.keyexit_rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) is False
    assert os.path.isdir(consensus_dir) is False
    # fp only
    run_file = os.path.join(results_dir, "run_consensus")
    with open(run_file, "w") as f:
        f.write("CONSAMBIG MOCK")
    G.keyexit_rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) is False
    assert os.path.isdir(consensus_dir) is False
    # DB extraction
    filename = "BLASTnontarget0.sequences"
    filepath = os.path.join(primer_qc_dir, filename)
    with open(filepath, "w") as f:
        f.write("DB extract Mock")
    assert os.path.isfile(filepath) is True
    G.keyexit_rollback("DB extraction", fp=filepath)
    assert os.path.isfile(filepath) is False
    # blast
    directory = blast_dir
    blastfile = "conserved_0_results.xml"
    filepath = os.path.join(blast_dir, blastfile)
    with open(filepath, "w") as f:
        f.write("BLAST MOCK")
    filename = blastfile
    assert os.path.isfile(filepath) is True
    G.keyexit_rollback("BLAST search", dp=directory, fn=filename)
    assert os.path.isfile(filepath) is False
    # primerblast
    directory = primerblast_dir
    primerblastfile = "primer_0_results.xml"
    filepath = os.path.join(primerblast_dir, primerblastfile)
    with open(filepath, "w") as f:
        f.write("PRIMERBLAST MOCK")
    filename = primerblastfile
    assert os.path.isfile(filepath) is True
    G.keyexit_rollback("BLAST search", dp=directory, fn=filename)
    assert os.path.isfile(filepath) is False
    # DB indexing
    dbfiles = [
        "Lb_curva.genomic", "Lb_curva.genomic.sqlite3.db",
        "Lb_curva.genomic.uni", "Lb_curva.genomic.2bit"]
    for files in dbfiles:
        filepath = os.path.join(primer_qc_dir, files)
        with open(filepath, "w") as f:
            f.write(files + " Mock")
        assert os.path.isfile(filepath) is True
    db_name = "Lb_curva.genomic"
    G.keyexit_rollback("DB indexing", dp=primer_qc_dir, search=db_name)
    for files in dbfiles:
        filepath = os.path.join(primer_qc_dir, files)
        assert os.path.isfile(filepath) is False
    # pangenome
    assert os.path.isdir(pangenome_dir) is True
    G.keyexit_rollback("pan-genome analysis", dp=pangenome_dir)
    assert os.path.isdir(pangenome_dir) is False
    shutil.rmtree(testdir)


if __name__ == "__main__":
    print(msg)
