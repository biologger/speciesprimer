#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import logging
import subprocess
import os
import csv
import concurrent.futures
from Bio import Entrez

pipe_bin = os.path.abspath(__file__)
pipe_dir = pipe_bin.split("pipeline")[0]
dict_path = os.path.join(pipe_dir, "dictionaries")

class GeneralFunctions:
    @staticmethod
    def logger(string_to_log):
        logging.info(time.strftime(
            "%d %b %Y %H:%M:%S: ", time.localtime())
            + str(string_to_log).strip())

    @staticmethod
    def run_subprocess(
            cmd, printcmd=True, logcmd=True, log=True, printoption=True):
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
                    if log:
                        logging.log(logging.INFO, output)
                    if printoption:
                        print(output)
                    elif printoption == "static":
                        print('\r ' + output, end='')
                else:
                    break

        while process.poll() is None:
            check_output()

    @staticmethod
    def run_shell(cmd, printcmd=True, logcmd=True, log=True, printoption=True):
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
                    if printoption:
                        print(output)
                    elif printoption == "static":
                        print('\r ' + output, end='')
                else:
                    break

        while process.poll() is None:
            check_output()

    @staticmethod
    def read_shelloutput(cmd, printcmd=False, logcmd=True, printoption=False):
        if logcmd:
            GeneralFunctions().logger("Run " + cmd)
        if printcmd:
            print("Run " + cmd)
        outputlist = []
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    outputlist.append(output)
                    if printoption:
                        print(output)
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
                writer.writerow(header)
            writer.writerows(inputlist)


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
            splitter = target.split("_")
            genus = splitter[0]
            species = splitter[1]
            sub = splitter[3][0:5]
            spec = species[0:5]
            try:
                geni = abbrev[genus]
            except KeyError:
                geni = genus[0:5]
            name = geni+"_"+spec+"_"+sub
        else:
            splitter = target.split("_")
            genus = splitter[0]
            species = splitter[1]
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
                info = "No taxid was found on NCBI"
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
