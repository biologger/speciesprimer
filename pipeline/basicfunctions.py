#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import logging
import subprocess
import os
import csv
import json
import shutil
import concurrent.futures
from Bio import Entrez
import tempfile
from collections import Counter
from datetime import timedelta

pipe_dir = os.path.dirname(os.path.abspath(__file__))
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')


class GeneralFunctions:
    @staticmethod
    def logger(string_to_log):
        logging.info(time.strftime(
            "%d %b %Y %H:%M:%S: ", time.localtime())
            + str(string_to_log).strip())

    @staticmethod
    def run_subprocess(
            cmd, printcmd=True, logcmd=True, printoption=True):
        if logcmd:
            GeneralFunctions().logger("Run " + " ".join(cmd))
        if printcmd:
            print("Run " + " ".join(cmd))
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    if printoption:
                        print(output)
                else:
                    break

        while process.poll() is None:
            check_output()

    @staticmethod
    def run_shell(cmd, printcmd=True, logcmd=True, log=True):
        if logcmd:
            GeneralFunctions().logger("Run " + cmd)
        if printcmd:
            print("\nRun " + cmd)
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    if log:
                        logging.log(logging.INFO, output)
                else:
                    break

        while process.poll() is None:
            check_output()

    @staticmethod
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

    @staticmethod
    def run_parallel(function, input_list, args=False, verbosity="bar"):
        outputlist = []
        total = len(input_list)
        start = None
        ProcessPool = concurrent.futures.ProcessPoolExecutor()
        bar = 0
        if verbosity == "bar":
            ending = 50 * " "
            print('\rprogress ' + str(0) + " % [" + str(ending) + "]", end='')
        with ProcessPool as executor:
            if args is False:
                future_seq = {
                    executor.submit(function, item):
                    item for item in input_list}
            else:
                future_seq = {
                    executor.submit(function, item, args):
                    item for item in input_list}

            for future in concurrent.futures.as_completed(future_seq):
                item = future_seq[future]
                try:
                    output = future.result()
                    outputlist.append(output)
                    if verbosity == "bar":
                        percent = round(100 / total * len(outputlist), 0)
                        if percent == percent // 1:
                            status = int(percent)
                            if status % 2 == 0:
                                bar = int(status / 2)
                                start = bar * "*"
                                ending = (50 - bar) * " "
                            print(
                                '\rprogress ' + str(status) + " % ["
                                + str(start) + str(ending) + "]", end='')
                except Exception as exc:
                    msg = (
                        '%r generated an exception: %s' % (item, exc)
                    )
                    print(msg)
                    GeneralFunctions().logger(msg)

        print("\n")
        return outputlist

    @staticmethod
    def create_directory(path_to_dir):
        if not os.path.isdir(path_to_dir):
            try:
                os.makedirs(path_to_dir)
            except OSError:
                if not os.path.isdir(path_to_dir):
                    raise

    @staticmethod
    def csv_writer(filepath, inputlist, header=None):
        with open(filepath, "w") as f:
            writer = csv.writer(f)
            if header:
                if any(isinstance(item, list) for item in header):
                    writer.writerows(header)
                else:
                    writer.writerow(header)
            writer.writerows(inputlist)

    @staticmethod
    def keyexit_rollback(stage, fp=None, dp=None, fn=None, search=None):
        logging.error("KeyboardInterrupt during " + stage, exc_info=True)
        msg = "No files affected"
        if dp and fn:
            fp = os.path.join(dp, fn)
            if os.path.isfile(fp):
                os.remove(fp)
                msg = "Remove " + fp
        elif dp and fp:
            if os.path.isdir(dp):
                shutil.rmtree(dp)
            else:
                dp = ""
            if os.path.isfile(fp):
                os.remove(fp)
            else:
                fp = ""
            msg = "Remove " + dp + " " + fp
        elif dp and search:
            deleted = []
            for files in os.listdir(dp):
                if files.startswith(search):
                    filepath = os.path.join(dp, files)
                    os.remove(filepath)
                    deleted.append(files)
            msg = "Remove " + " ".join(deleted)
        elif fp:
            if os.path.isfile(fp):
                os.remove(fp)
        else:
            if os.path.isdir(dp):
                shutil.rmtree(dp)
                msg = "Remove " + dp
        print("\n" + msg + "\n")
        GeneralFunctions().logger(msg)


class HelperFunctions:

    @staticmethod
    def subspecies_handler(target, mode="underscore"):
        GeneralFunctions().logger(
            "Run: subspecies_handler(" + target + " , " + mode + ")")
        if "_subsp_" in target:
            if mode == "underscore":
                species = (
                    target.split("_")[1] + "_subsp_"
                    + target.split("_")[3])
            if mode == "space":
                species = (
                    target.split("_")[1] + " subsp. "
                    + target.split("_")[3])
        else:
            species = target.split("_")[1]
        return species

    @staticmethod
    def create_non_target_list(target):
        nontarget_list = []
        genus = target.split("_")[0]
        species = HelperFunctions().subspecies_handler(target, "space")
        spec_list = os.path.join(dict_path, "species_list.txt")
        with open(spec_list, "r") as species_list:
            for line in species_list.readlines():
                line = line.strip()
                if 'subs' in species:
                    if not genus+" "+species == str(line):
                        nontarget_list.append(line)
                else:
                    if genus+" "+species not in str(line):
                        nontarget_list.append(line)

            return nontarget_list

    @staticmethod
    def abbrev(target):
        abbrev = {}
        with open(os.path.join(dict_path, "genus_abbrev.csv")) as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            for row in reader:
                species = row[0]
                short = row[1]
                abbrev.update({species: short})
        if "subsp" in target:
            genus = target.split("_")[0]
            species = target.split("_")[1]
            sub = target.split("_")[3][0:5]
            spec = species[0:5]
            try:
                geni = abbrev[genus]
            except KeyError:
                geni = genus[0:5]
            name = geni+"_"+spec+"_"+sub
        else:
            genus = target.split("_")[0]
            species = target.split("_")[1]
            spec = species[0:5]
            try:
                geni = abbrev[genus]
            except KeyError:
                geni = genus[0:5]
            name = geni+"_"+spec
        return name

    @staticmethod
    def check_input(target, email):
        GeneralFunctions().logger("Run: check_input()")
        try:
            Entrez.email = email
            searchtaxid = Entrez.esearch(db="taxonomy", term=target)
            taxidresult = Entrez.read(searchtaxid)
            taxid = taxidresult["IdList"]
            if len(taxid) == 1:
                return taxid[0]
            else:
                error = taxidresult['ErrorList']
                info = "No taxid was found on NCBI\nError: " + str(error)
                print(info)
                GeneralFunctions().logger("> " + info)
        except OSError:
            info = (
                "ERROR: Taxid for " + target
                + " not found, please check spelling and internet connection")
            print(info)
            GeneralFunctions().logger("> " + info)
            time.sleep(2)
            raise

    @staticmethod
    def check_species_syn(taxid, email, target):
        Entrez.email = email
        searchsyn = Entrez.efetch(db="taxonomy", id=taxid)
        synresult = Entrez.read(searchsyn)
        scienctificname = synresult[0]['ScientificName']
        synonym = synresult[0]['OtherNames']['Synonym']
        includes = synresult[0]['OtherNames']['Includes']
        synonyms = synonym + includes
        if synonyms == []:
            return None
        else:
            synwarn = []
            target_name = " ".join(target.split("_"))
            if not target_name == scienctificname:
                synwarn.append(scienctificname)
            for item in synonyms:
                if not item == target_name:
                    synwarn.append(item)
            if synwarn == []:
                return None
            else:
                info = ("Warning synonyms for this species were found...")
                info2 = ("Adding synonyms to exception in config.json.")
                print("\n" + info)
                print(synwarn)
                print(info2 + "\n")
                GeneralFunctions().logger("> " + info)
                GeneralFunctions().logger(synwarn)
                GeneralFunctions().logger("> " + info2)
                return synwarn

    @staticmethod
    def get_email_for_Entrez(email=None):
        if os.path.isfile(tmp_db_path):
            with open(tmp_db_path) as f:
                for line in f:
                    tmp_db = json.loads(line)
            if email:
                if "@" in email and "." in email:
                    tmp_db.update({'email': email})
                    with open(tmp_db_path, 'w') as f:
                        f.write(json.dumps(tmp_db))
            try:
                mail = tmp_db['email']
                if "@" in mail and "." in mail:
                    email = mail.strip()
            except KeyError:
                email = input(
                    "To make use of NCBI's E-utilities, "
                    "Please enter your email address. \n")
                if "@" in email and "." in email:
                    with open(tmp_db_path) as f:
                        for line in f:
                            tmp_db = json.loads(line)
                    tmp_db.update({'email': email})
                    with open(tmp_db_path, 'w') as f:
                        f.write(json.dumps(tmp_db))
                else:
                    print("Not a valid email adress")
                    HelperFunctions().get_email_for_Entrez()
        else:
            if email:
                pass
            else:
                email = input(
                    "To make use of NCBI's E-utilities, "
                    "Please enter your email address. \n")

            if "@" in email and "." in email:
                tmp_db = {
                    'email': email, 'new_run': {
                        'path': '', 'targets': {},
                        'same_settings': False,
                        'change_settings': False
                    }
                }
                with open(tmp_db_path, 'w') as f:
                    f.write(json.dumps(tmp_db))
            else:
                print("Not a valid email adress")
                HelperFunctions().get_email_for_Entrez()

        return email


class ParallelFunctions:

    @staticmethod
    def get_seq_fromDB(extractdata, db):
        [accession, start, stop] = extractdata
        fasta = []
        seq_cmd = " ".join([
            "blastdbcmd", "-db", db, "-entry", str(accession),
            "-range", str(start) + "-" + str(stop), "-outfmt", "%f"])
        while fasta == []:
            fasta = GeneralFunctions().read_shelloutput(seq_cmd)
        return fasta

    @staticmethod
    def MFEprimer_template(primerinfo, args):
        [primer_qc_dir, mfethreshold] = args
        result = []
        [nameF, seqF, nameR, seqR, templ_seq] = primerinfo
        with tempfile.NamedTemporaryFile(
            mode='w+', dir=primer_qc_dir, prefix="primer",
            suffix=".fa", delete=False
        ) as primefile:
            primefile.write(
                ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")

        db = "template.sequences"
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + db
            + " -k 9 --tab --ppc 10")
        while result == []:
            result = GeneralFunctions().read_shelloutput(cmd)
        os.unlink(primefile.name)
        if len(result) == 2:
            val = result[1].split("\t")
            pp_F = "_".join(val[1].split("_")[0:-1])
            pp_R = "_".join(val[2].split("_")[0:-1])
            p_F = "_".join(val[1].split("_")[0:-2])
            primername = val[3]
            ppc = float(val[4])
            if (
                pp_F == pp_R and p_F == primername
                and ppc >= float(mfethreshold)
            ):
                ppc_val = ppc - float(mfethreshold)
                return [[nameF, seqF, nameR, seqR, templ_seq, ppc_val], result]
            else:
                return [[None], result]
        else:
            return [[None], result]

    @staticmethod
    def MFEprimer_nontarget(primerinfo, args):
        result = []
        nameF, seqF, nameR, seqR, templ_seq, ppc_val = primerinfo
        [dbfilepath, primer_qc_dir] = args
        dbfile = os.path.basename(dbfilepath)
        with tempfile.NamedTemporaryFile(
            mode='w+', dir=primer_qc_dir, prefix="primer",
            suffix=".fa", delete=False
        ) as primefile:
            primefile.write(
                ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + dbfile
            + " -k 9 --tab --ppc 10")
        while result == []:
            result = GeneralFunctions().read_shelloutput(cmd)
        os.unlink(primefile.name)
        if not len(result) == 1:
            for index, item in enumerate(result):
                if index > 0:
                    val = item.split("\t")
                    result_ppc = float(val[4])
                    if result_ppc > ppc_val:
                        return [[None], result]

        return [primerinfo, result]

    @staticmethod
    def MFEprimer_assembly(primerinfo, args):
        [primer_qc_dir, db, mfethreshold] = args
        result = []
        target_product = []
        nameF, seqF, nameR, seqR, templ_seq, ppc_val = primerinfo
        with tempfile.NamedTemporaryFile(
            mode='w+', dir=primer_qc_dir, prefix="primer",
            suffix=".fa", delete=False
        ) as primefile:
            primefile.write(
                ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + db
            + " -k 9 --tab --ppc 10")
        while result == []:
            result = GeneralFunctions().read_shelloutput(cmd)
        os.unlink(primefile.name)
        for index, item in enumerate(result):
            if index > 0:
                val = item.split("\t")
                result_ppc = float(val[4])
                product_len = int(val[5])
                targetID = val[3]
                if result_ppc == ppc_val + mfethreshold:
                    target_product.append(targetID)
                elif result_ppc > ppc_val:
                    return [[None], result]
        counts = Counter(target_product)
        for item in counts.keys():
            if counts[item] == 1:
                return [primerinfo, result]
            else:
                return [[None], result]

    @staticmethod
    def index_database(inputfilepath):
        primer_qc_dir = os.path.dirname(inputfilepath)
        db_name = os.path.basename(inputfilepath)
        db_path = inputfilepath + ".sqlite3.db"
        if os.path.isfile(db_path) is True:
            msg = " ".join([db_name, "DB already exists"])
            print(msg)
            GeneralFunctions().logger(msg)
            return 0
        elif os.stat(inputfilepath).st_size == 0:
            db_name = os.path.basename(inputfilepath)
            msg = " ".join(["Problem with", db_name, "input file is empty"])
            GeneralFunctions().logger("> " + msg)
            print("\n!!!" + msg + "!!!\n")
            os.remove(inputfilepath)
            return msg
        else:
            GeneralFunctions().logger("> Start index non-target DB " + db_name)
            print("\nStart index " + db_name)
            start = time.time()
            cmd = "IndexDb.sh " + inputfilepath + " 9"
            try:
                GeneralFunctions().run_shell(
                        cmd, printcmd=True, logcmd=True, log=False)
            except (KeyboardInterrupt, SystemExit):
                GeneralFunctions().keyexit_rollback(
                        "DB indexing", dp=primer_qc_dir, search=db_name)
                raise
            end = time.time() - start
            GeneralFunctions().logger(
                "Run: index_Database(" + db_name + ") time: "
                + str(timedelta(seconds=end)))
            print("Done indexing " + db_name)
            GeneralFunctions().logger("> Done indexing " + db_name)
            return 0
