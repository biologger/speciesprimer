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
