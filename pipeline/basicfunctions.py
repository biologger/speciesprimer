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
    def comm_log(string_to_log, newline=False, output=None):
        GeneralFunctions().logger(string_to_log)
        if output:
            with output:
                if newline:
                    print("\n" + string_to_log + "\n")
                else:
                    print(string_to_log)
        else:
            if newline:
                print("\n" + string_to_log + "\n")
            else:
                print(string_to_log)




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
            stderr=subprocess.STDOUT)

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
    def accession_from_filename(filename, version=True):
        if "GCF" in filename or "GCA" in filename:
            accession = "_".join(filename.split("_")[0:2])
            if version:
                accession = "v".join(accession.split("."))
        else:
            accession = "_".join(filename.split("_")[0:-1])
        return accession

    @staticmethod
    def genomicversion_from_accession(accession):
        if "GCF" in accession or "GCA" in accession:
            accession = "_".join(accession.split("_")[0:2])
            accession = ".".join(accession.split("v"))
        return accession

    @staticmethod
    def advanced_pipe_config(path_to_configfile):
        options = [
            ["genus_abbrev", os.path.join(dict_path, "genus_abbrev.csv")],
            ["species_list", os.path.join(dict_path, "species_list.txt")],
            ["p3settings", os.path.join(dict_path, "p3parameters")],
            ["excludedgis", os.path.join(dict_path, "no_blast.gi")]]

        with open(path_to_configfile) as f:
            for line in f:
                adv_sett = json.loads(line)

        if "certificate" in adv_sett.keys():
            certfile = adv_sett['certificate']
            if os.path.isfile(certfile):
                filename_ext = os.path.basename(certfile)
                filename, file_ext = os.path.splitext(filename_ext)
                if file_ext == ".crt" or file_ext == ".cer":
                    certloc = os.path.join(
                        "/", "usr", "local", "share", "ca-certificates",
                        filename + ".crt")
                    shutil.copy(certfile, certloc)
                    GeneralFunctions().run_subprocess(
                            ["update-ca-certificates"])
                    certcmd = [
                        "echo", "export",
                        "REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt",
                        ">>" ,"/root/.bashrc"]
                    GeneralFunctions().run_shell(" ".join(certcmd))
                    GeneralFunctions().run_shell(
                            " ".join(["source", "/root/.bashrc"]))
                else:
                    msg = "Certificates in .crt or .cer format are required"
                    print(msg)
                    GeneralFunctions().logger(" ".join(msg))
                    return 1
            else:
                msg = " ".join(["Certificate", certfile, "not found"])
                print(msg)
                GeneralFunctions().logger(" ".join(msg))
                return 1

        for option in options:
            if option[0] in adv_sett.keys():
                settings = adv_sett[option[0]]
                print(settings)
                if isinstance(settings, list):
                    if isinstance(settings[0], list):
                        GeneralFunctions().csv_writer(option[1], settings)
                    else:
                        with open(option[1], "w") as f:
                            for item in settings:
                                f.write(str(item) + "\n")
                else:
                    if os.path.isfile(settings):
                        filename, file_ext = os.path.splitext(option[1])
                        refname, ref_ext = os.path.splitext(settings)
                        if ref_ext == file_ext:
                            os.remove(option[1])
                            shutil.copy(settings, option[1])
                        else:
                            msg = [
                                "A", file_ext, "file is required as",
                                option[0], "settings file"]
                            print(" ".join(msg))
                            GeneralFunctions().logger(" ".join(msg))
                            return 1
                    else:
                        msg = (
                            "A list containing the settings or a complete "
                            "settings file is required")
                        print(msg)
                        GeneralFunctions().logger(msg)
                        return 1
        return 0


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
        elif "_pv_" in target:
            if mode == "underscore":
                species = (
                    target.split("_")[1] + "_pv_"
                    + target.split("_")[3])
            if mode == "space":
                species = (
                    target.split("_")[1] + " pv. "
                    + target.split("_")[3])
        else:
            if len(target.split("_")) > 1:
                species = target.split("_")[1]
            else:
                species = target
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
            if len(target.split("_")) == 1:
                return target
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
            taxidresult = Entrez.read(searchtaxid, validate=False)
            taxid = taxidresult["IdList"]
            if len(taxid) == 1:
                syn = HelperFunctions().check_species_syn(
                                                    taxid[0], email, target)
                return taxid[0], syn

            error = taxidresult['ErrorList']
            info = "No taxid was found on NCBI\nError: " + str(error)
            print(info)
            GeneralFunctions().logger("> " + info)
            return None, None
        except OSError:
            info = (
                "ERROR: Taxid for " + target
                + " not found, please check internet connection")
            print(info)
            GeneralFunctions().logger("> " + info)
            time.sleep(2)
            raise

    @staticmethod
    def check_species_syn(taxid, email, target):
        Entrez.email = email
        try:
            searchsyn = Entrez.efetch(db="taxonomy", id=taxid)
            synresult = Entrez.read(searchsyn, validate=False)
            scienctificname = synresult[0]['ScientificName']
            synonym = synresult[0]['OtherNames']['Synonym']
            includes = synresult[0]['OtherNames']['Includes']
            equivalents = synresult[0]['OtherNames']['EquivalentName']

            synonyms = synonym + includes + equivalents
            if synonyms != []:
                synwarn = []
                target_name = " ".join(target.split("_"))
                if not target_name == scienctificname:
                    synwarn.append(scienctificname)
                for item in synonyms:
                    if not item == target_name:
                        synwarn.append(item)
                if synwarn != []:
                    info = ("Warning synonyms for this species were found...")
                    info2 = ("Adding synonyms to exception in config.json.")
                    print("\n" + info)
                    print(synwarn)
                    print(info2 + "\n")
                    GeneralFunctions().logger("> " + info)
                    GeneralFunctions().logger(synwarn)
                    GeneralFunctions().logger("> " + info2)
                    return synwarn
            return None
        except OSError:
            info = (
                "SpeciesPrimer is unable to connect to the Entrez server,"
                " please check the internet connection and try again later")
            print(info)
            logging.error("> " + info, exc_info=True)
            from speciesprimer import errors
            errors.append([target, info])
            raise

    @staticmethod
    def get_email_for_Entrez(email=None):

        def tmp_db_email(email):
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
            return email

        def user_input_email(email):
            if not email:
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

        if os.path.isfile(tmp_db_path):
            email = tmp_db_email(email)
        else:
            email = user_input_email(email)

        return email


    @staticmethod
    def BLASTDB_check(config):
        if config.customdb:
            DBname = config.customdb
        else:
            DBname = "nt"
        if (
            os.path.isfile(DBname + ".nsq") is True or
            os.path.isfile(DBname + ".nal") is True or
            os.path.isfile(
                os.path.join("/", "blastdb", DBname + ".nal")) is True
        ):
            return 0

        msg = "No BLAST DB was found " + DBname
        raise BlastDBError(msg)


class BlastDBError(Exception):
    pass


class ParallelFunctions:

    @staticmethod
    def get_seq_fromDB(extractdata, db):
        [accession, start, stop] = extractdata
        fasta = []
        seq_cmd = [
            "blastdbcmd", "-db", db, "-entry", str(accession),
            "-range", str(start) + "-" + str(stop), "-outfmt", "%f"]
        while fasta == []:
            fasta = GeneralFunctions().read_shelloutput(seq_cmd)
        return fasta

    def MFEprimer_run(primerinfo, args):
        result = []
        name, seqF, seqR = primerinfo
        [dbfilepath, primer_qc_dir] = args
        with tempfile.NamedTemporaryFile(
            mode='w+', dir=primer_qc_dir, prefix="primer",
            suffix=".fa", delete=False
        ) as primefile:
            primefile.write(
                ">" + name + "_F\n" + seqF + "\n>" + name + "_R\n" + seqR + "\n")
        cmd = [
            "MFEprimer.py", "-i", primefile.name, "-d", dbfilepath,
            "-k", "9", "--tab", "--ppc", "10"]

        while result == []:
            result = GeneralFunctions().read_shelloutput(cmd)
        os.unlink(primefile.name)
        return result


    @staticmethod
    def MFEprimer_singleton(primerinfo, args):
        [primer_qc_dir, db, mfethreshold, short] = args
        nameF, seqF, seqR, ppc_val = primerinfo
        targetname = "_".join(nameF.split(short)[1].split("_")[0:-3])
        dbname = os.path.basename(os.path.dirname(db))
        if dbname == targetname:
            result = ParallelFunctions().MFEprimer_assembly(
                                                        primerinfo, args[0:3])
        else:
            result = ParallelFunctions().MFEprimer_nontarget(
                    primerinfo, [db, primer_qc_dir])
        return result


    @staticmethod
    def index_database(inputfilepath):
        primer_qc_dir = os.path.dirname(inputfilepath)
        db_name = os.path.basename(inputfilepath)
        db_path = inputfilepath + ".sqlite3.db"
        if os.path.isfile(db_path) is True:
            msg = " ".join([db_name, "DB already exists"])
            GeneralFunctions().comm_log(msg, True)
            if os.stat(db_path).st_size == 0:
                msg = " ".join(["Problem with", db_name, "db file is empty"])
                GeneralFunctions().comm_log("> " + msg)
                os.remove(db_path)
            else:
                return 0

        if os.stat(inputfilepath).st_size == 0:
            db_name = os.path.basename(inputfilepath)
            msg = " ".join(["Problem with", db_name, "input file is empty"])
            GeneralFunctions().comm_log("> " + msg, True)
            os.remove(inputfilepath)
            return msg


        GeneralFunctions().logger("> Start index non-target DB " + db_name)
        print("\nStart index " + db_name)
        start = time.time()
        cmd = ["IndexDB.py", inputfilepath, "-k", "9"]
        try:
            GeneralFunctions().run_subprocess(
                    cmd, True, True, False)
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
