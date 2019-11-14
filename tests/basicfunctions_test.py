#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import pytest
import json

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
sys.path.append(pipe_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
dbpath = os.path.join(tmpdir, "customdb.fas")

#            import traceback
#            traceback.print_exc()

from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G

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
            self.repeat +=1
        elif self.repeat < 2:
            prompt_dict = {user_email: "hfkjdsa@lgfhlghcom"}
            assert read_email_from_db() == self.outcome
            self.repeat +=1
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
    else:
        return inputitem

def test_run_parallel_exceptions():
    inputlist = [
            "lactis_subsp_lactis", "lactis subsp. lactis",
            "enterica", "enterica"]
    result = G.run_parallel(a_function, inputlist, args=False, verbosity="")
    assert result == ['lactis_subsp_lactis', 'enterica', 'enterica']

def test_subspecies_handler():
    outcome = [
            "lactis_subsp_lactis", "lactis subsp. lactis",
            "enterica", "enterica"]
    targets = ["Lactococcus_lactis_subsp_lactis", "Salmonella_enterica"]
    modes = ["underscore", "space"]
    outnum = 0
    for target in targets:
        for mode in modes:
            species = H.subspecies_handler(target, mode=mode)
            assert species == outcome[outnum]
            outnum += 1

def test_check_input_fail(monkeypatch):
    taxid = H.check_input("Lactobacius_curvatus", "biologger@protonmail.com")
    assert taxid == None

def test_create_non_target_list():
    targetdef = []
    target = "Lactococcus_lactis_subsp_lactis"
    nontargetlist = H.create_non_target_list(target)
    for item in nontargetlist:
        if (
            item == "Lactococcus lactis subsp. lactis" or
            item == "Lactococcus lactis subsp lactis"):
            targetdef.append(item)
    assert targetdef == []

def test_rollback():
    testdir = os.path.join("/", "primerdesign", "tests")
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
    assert os.path.isdir(os.path.join(target_dir, annotdir)) == True
    G.rollback("annotation", dp=os.path.join(target_dir, outdir))
    assert os.path.isdir(os.path.join(target_dir, annotdir)) == False
    # p3
    p3_file = os.path.join(primer_dir, "primer3_output")
    output_file = p3_file
    with open(p3_file, "w") as f:
        f.write("PRIMER3MOCK")
    assert os.path.isfile(p3_file) == True
    G.rollback("primer3 run", fp=output_file)
    assert os.path.isfile(p3_file) == False
    # alignments
    run_file = os.path.join(results_dir, "run_prank")
    with open(run_file, "w") as f:
        f.write("PRANK MSA MOCK")
    assert os.path.isfile(run_file) == True
    G.rollback("Prank MSA run", dp=alignments_dir, fp=run_file)
    assert os.path.isfile(run_file) == False
    assert os.path.isdir(alignments_dir) == False
    # consensus
    run_file = os.path.join(results_dir, "run_consensus")
    with open(run_file, "w") as f:
        f.write("CONSAMBIG MOCK")
    assert os.path.isfile(run_file) == True
    G.rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) == False
    assert os.path.isdir(consensus_dir) == False
    # dp only
    run_file = os.path.join(results_dir, "run_consensus")
    G.create_directory(consensus_dir)
    assert os.path.isfile(run_file) == False
    G.rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) == False
    assert os.path.isdir(consensus_dir) == False
    # fp only
    run_file = os.path.join(results_dir, "run_consensus")
    with open(run_file, "w") as f:
        f.write("CONSAMBIG MOCK")
    G.rollback("consensus run", dp=consensus_dir, fp=run_file)
    assert os.path.isfile(run_file) == False
    assert os.path.isdir(consensus_dir) == False
    # DB extraction
    filename = "BLASTnontarget0.sequences"
    filepath = os.path.join(primer_qc_dir, filename)
    with open(filepath, "w") as f:
        f.write("DB extract Mock")
    assert os.path.isfile(filepath) == True
    G.rollback("DB extraction", fp=filepath)
    assert os.path.isfile(filepath) == False
    # blast
    directory = blast_dir
    blastfile = "conserved_0_results.xml"
    filepath = os.path.join(blast_dir, blastfile)
    with open(filepath, "w") as f:
        f.write("BLAST MOCK")
    filename = blastfile
    assert os.path.isfile(filepath) == True
    G.rollback("BLAST search", dp=directory, fn=filename)
    assert os.path.isfile(filepath) == False
    # primerblast
    directory = primerblast_dir
    primerblastfile = "primer_0_results.xml"
    filepath = os.path.join(primerblast_dir, primerblastfile)
    with open(filepath, "w") as f:
        f.write("PRIMERBLAST MOCK")
    filename = primerblastfile
    assert os.path.isfile(filepath) == True
    G.rollback("BLAST search", dp=directory, fn=filename)
    assert os.path.isfile(filepath) == False
    # DB indexing
    dbfiles = [
        "Lb_curva.genomic", "Lb_curva.genomic.sqlite3.db",
        "Lb_curva.genomic.uni", "Lb_curva.genomic.2bit"]
    for files in dbfiles:
        filepath = os.path.join(primer_qc_dir, files)
        with open(filepath, "w") as f:
            f.write(files + " Mock")
        assert os.path.isfile(filepath) == True
    db_name = "Lb_curva.genomic"
    G.rollback("DB indexing", dp=primer_qc_dir, search=db_name)
    for files in dbfiles:
        filepath = os.path.join(primer_qc_dir, files)
        assert os.path.isfile(filepath) == False
    # pangenome
    assert os.path.isdir(pangenome_dir) == True
    G.rollback("pan-genome analysis", dp=pangenome_dir)
    assert os.path.isdir(pangenome_dir) == False

    shutil.rmtree(testdir)
