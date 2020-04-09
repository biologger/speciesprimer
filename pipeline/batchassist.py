#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import sys
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

# paths
pipe_dir = os.path.dirname(os.path.abspath(__file__))
dict_path = os.path.join(pipe_dir, "dictionaries")


class Input:
    def __init__(self):
        self.target_list = []
        self.config_dict = {}
        self.input_dict = {}
        self.config_paths = []

    def gui_runner(self, mode, data):
        if mode == "continue":
            path, targets = data
            config_dict = Output().run_gui_output(targets, path)
            return config_dict

        if mode == "new":
            self.config_dict = data
            for target in self.config_dict.keys():
                self.target_list.append(target)
                self.write_config_file(target)
            return self.config_dict

    def initiate(self):
        mode = input(
            "Create new config files or start pipeline with previously "
            "generated files?\ntype (n)ew or (s)tart:\n")
        if (mode.lower() == "n" or mode.lower() == "new"):
            config_dict = self.main()
            return config_dict
        if (mode.lower() == "s" or mode.lower() == "start"):
            config_dict = Output().run_output()
            return config_dict

        print("no valid input, type: (n)ew or (s)tart")
        sys.exit()

    def targetinput(self):
        targets = input(
            "Please specify target species (comma separated) "
            "or type help for format examples: \n")
        return targets

    def helpmessage(self, targets):
        print(
            "Accepted input formats are e.g. Lactobacillus_curvatus, "
            "Lactobacillus delbrueckii, Lactococcus lactis subsp lactis, "
            "Lactococcus_lactis_subsp_cremoris")
        targets = input(
            "Please specify target species (comma separated) "
            "or type help for format examples: \n")
        return targets

    def parse_targets(self, targets):
        if isinstance(targets, list):
            items = targets
        else:
            items = list(targets.split(","))
        for i, item in enumerate(items):
            x = item.strip()
            x = x.capitalize()
            self.target_list.insert(i, x)

        for i, item in enumerate(self.target_list):
            target = item.strip()
            if " " in target:
                self.target_list.remove(target)
                new_target = "_".join(target.split(" "))
                self.target_list.insert(i, new_target)
        for target in self.target_list:
            self.config_dict.update({target: {"target": target}})
        return self.target_list

    def get_path(self, target, index, listlen):
        if "path" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            inpath = input(
                "Please specify a path to use as the working directory "
                "or hit return to use the current working directory:\n")
            if inpath:
                if os.path.isabs(inpath):
                    path = inpath
                else:
                    path = os.path.join(os.getcwd(), inpath)
            else:
                path = os.getcwd()
            if index == 0:
                if not self.value_for_all("path", path, listlen):
                    self.config_dict[target].update({"path": path})
            else:
                self.config_dict[target].update({"path": path})

            print("path", path)

    def work_offline(self, target, index, listlen):
        if "offline" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            limit = input(
                "Work offline with local genome assemblies?"
                "\n(y)es/(n)o, default=(n)\n> ")
            if limit.lower() == "y" or limit.lower() == "yes":
                offline = True
                skip_download = True
                assembly_list = ["offline"]

            else:
                offline = False
                skip_download = False
            if index == 0:
                if not self.value_for_all("offline", offline, listlen):
                    self.config_dict[target].update({"offline": offline})
                    if skip_download:
                        self.config_dict[target].update(
                            {"skip_download": skip_download})
                        self.config_dict[target].update(
                            {"assemblylevel": assembly_list})

                else:
                    if skip_download:
                        for selectedtarget in self.target_list:
                            self.config_dict[selectedtarget].update(
                                {"skip_download": skip_download})
                            self.config_dict[selectedtarget].update(
                                {"assemblylevel": assembly_list})

            else:
                self.config_dict[target].update({"offline": offline})
                if skip_download:
                    self.config_dict[target].update(
                        {"skip_download": skip_download})
                    self.config_dict[target].update(
                        {"assemblylevel": assembly_list})

            print("offline", offline)

    def get_skip_download(self, target, index, listlen):
        if "skip_download" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            skip_download = input(
                "Skip the download of Genomes from NCBI?"
                "\n(y)es/(n)o, default=(n)\n> ")
            if skip_download.lower() == "y" or skip_download.lower() == "yes":
                skip_download = True
                assembly_list = ["offline"]
            else:
                skip_download = False

            if index == 0:
                if not self.value_for_all(
                    "skip_download", skip_download, listlen
                ):
                    self.config_dict[target].update(
                        {"skip_download": skip_download})
                    if skip_download:
                        self.config_dict[target].update(
                            {"assemblylevel": assembly_list})
                else:
                    if skip_download:
                        for selectedtarget in self.target_list:
                            self.config_dict[selectedtarget].update(
                                {"assemblylevel": assembly_list})
            else:
                self.config_dict[target].update(
                    {"skip_download": skip_download})
                if skip_download:
                    self.config_dict[target].update(
                        {"assemblylevel": assembly_list})
            print("skip_download", skip_download)

    def get_customdb(self, target, index, listlen):
        if "customdb" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            customdb = input(
                "Do you want to use a custom database for blastn?\n"
                "Specifiy the absolute filepath of the custom database "
                "(e.g. '/blastdb/nontarget.fasta' "
                "or '/blastdb/ref_prok_rep_genomes') "
                "or hit return to skip, default=None\n> ")
            if isinstance(customdb, str) and len(customdb) > 0:
                if (
                    os.path.isfile(customdb + ".nsq") is True or
                    os.path.isfile(customdb + ".nal") is True or
                    os.path.isfile(
                        os.path.join(
                                "/", "blastdb", customdb + ".nal")) is True
                ):
                    customdb = customdb
                else:
                    print("No BLAST DB was found " + customdb)
                    return self.get_customdb(target, index, listlen)
            else:
                customdb = None

            if index == 0:
                if not self.value_for_all("customdb", customdb, listlen):
                    self.config_dict[target].update(
                        {"customdb": customdb})
            else:
                self.config_dict[target].update({"customdb": customdb})
            print("customdb", customdb)

    def get_exception(self, target, index, listlen):
        if "exception" not in self.config_dict[target].keys():
            exception = []
            print("\n" + target + ":")
            inexception = input(
                "Primer binding to this non-target species is tolerated.\n"
                "Provide a species name or hit return to skip:\n> ")
            if inexception:
                inputexcept = list(inexception.split(","))
                for i, item in enumerate(inputexcept):
                    x = item.strip()
                    x = x.capitalize()
                    exception.insert(i, x)
                if index == 0:
                    if not self.value_for_all("exception", exception, listlen):
                        self.config_dict[target].update(
                                                    {"exception": exception})
                else:
                    self.config_dict[target].update({"exception": exception})
                print("exception", exception)
            else:
                if index == 0:
                    if not self.value_for_all("exception", exception, listlen):
                        self.config_dict[target].update(
                                                    {"exception": exception})
                else:
                    self.config_dict[target].update({"exception": exception})
                print("exception", exception)

    def get_single(self, target, index, listlen):
        if "strains" not in self.config_dict[target].keys():
            strains = []
            if "strain" in self.config_dict[target]["runmode"]:
                print("\n" + target + ":")
                insingle = input(
                    "Start of filename(s) of annotated fna file, e.g. "
                    "GCF_XYZXYZXYZv1, will only search for singletons "
                    "(strain specific primers) for "
                    "specified genome(s). (comma separated), default=[]")

                inputsingle = list(insingle.split(","))
                for i, item in enumerate(inputsingle):
                    x = item.strip()
                    strains.insert(i, x)

                self.config_dict[target].update({"strains": strains})
                print("strains", strains)
            else:
                self.config_dict[target].update({"strains": strains})

    def value_for_all(self, key, value, listlen):
        if listlen > 1:
            forall = input(
                "Use this value for all targets?\n(y)es/(n)o, default=(y)\n> ")
            if forall.lower() == "n" or forall.lower() == "no":
                return False

            for target in self.target_list:
                self.config_dict[target].update({key: value})
            return True

        return False

    def write_config_file(self, target):
        path = self.config_dict[target]["path"]
        target_dir = os.path.join(path, target)
        genomic_dir = os.path.join(target_dir, "genomic_fna")
        config_dir = os.path.join(target_dir, "config")
        config_file = os.path.join(config_dir, "config.json")
        G.create_directory(target_dir)
        G.create_directory(config_dir)
        G.create_directory(genomic_dir)
        self.config_paths.append(config_file)
        with open(config_file, "w") as f:
            f.write(json.dumps(self.config_dict[target]))

    def boolean_evaluation(self, userinput):
        return userinput.lower() == "y" or userinput.lower() == "yes"

    def number_evaluation(self, userinput, settings):
        [target, index, listlen, key] = settings
        evalnum = self.input_dict[key]["eval"][1][0]
        if evalnum == "int":
            ntype = int
        else:
            ntype = float
        if userinput:
            try:
                value = ntype(userinput)
            except ValueError:
                print("No valid input of type " + evalnum)
                return self.evaluation(target, index, listlen, key)
        else:
            value = self.input_dict[key]["eval"][1][1]
        return value

    def option_evaluation(self, userinput, settings):
        [target, index, listlen, key] = settings
        options = self.input_dict[key]["eval"][1]
        if userinput:
            try:
                valid = int(userinput)
            except ValueError:
                print("No valid input of type", int)
                return self.evaluation(target, index, listlen, key)

            if valid in options:
                value = valid
            else:
                print("No valid input, see options", options)
                return self.evaluation(target, index, listlen, key)
        else:
            value = self.input_dict[key]["eval"][1][0]
        return value

    def options_evaluation(self, userinput, settings):
        [target, index, listlen, key] = settings
        strippedlist = []
        options = self.input_dict[key]["eval"][1]
        if userinput:
            items = list(userinput.split(","))
            for i, item in enumerate(items):
                x = item.strip()
                strippedlist.insert(i, x)
            value = []
            for s_item in strippedlist:
                if s_item in options:
                    value.append(s_item)
                else:
                    print("\nNo valid input:")
                    print(s_item)
                    return self.evaluation(target, index, listlen, key)
        else:
            value = [self.input_dict[key]["eval"][1][0]]
        return value

    def evaluation(self, target, index, listlen, key):
        settings = [target, index, listlen, key]
        userinput = input(self.input_dict[key]["prompt"])
        evaltype = self.input_dict[key]["eval"][0]
        if evaltype == "boolean":
            value = self.boolean_evaluation(userinput)
        if evaltype == "number":
            value = self.number_evaluation(userinput, settings)
        if evaltype == "option":
            value = self.option_evaluation(userinput, settings)
        if evaltype == "options":
            value = self.options_evaluation(userinput, settings)
        return value

    def get_userinput(self, target, index, listlen, key):
        if key not in self.config_dict[target].keys():
            print("\n" + target + ":")
            newvalue = self.evaluation(target, index, listlen, key)
            if index == 0:
                if not self.value_for_all(key, newvalue, listlen):
                    self.config_dict[target].update({key: newvalue})
            else:
                self.config_dict[target].update({key: newvalue})
            print(key, newvalue)

    def main(self):
        targets = self.targetinput()
        while (targets == "help" or targets == "" or targets is None):
            targets = self.helpmessage(targets)

        dictpath = os.path.join(
                dict_path, "default", "batchassist_inputdict.json")
        with open(dictpath) as f:
            for line in f:
                self.input_dict = json.loads(line)
        setlist = [
            "runmode", "blastseqs", "virus", "genbank",
            "assemblylevel", "qc_gene", "ignore_qc",
            "skip_tree", "minsize", "maxsize", "probe", "mfold", "mpprimer",
            "mfethreshold", "evalue", "nuc_identity", "intermediate", "nolist"]

        self.parse_targets(targets)
        listlen = len(self.target_list)
        # get input
        for i, target in enumerate(self.target_list):
            self.get_path(target, i, listlen)
            self.work_offline(target, i, listlen)
            self.get_skip_download(target, i, listlen)
            self.get_exception(target, i, listlen)
            self.get_customdb(target, i, listlen)
            for item in setlist:
                self.get_userinput(target, i, listlen, item)
            self.get_single(target, i, listlen)

        for target in self.target_list:
            self.write_config_file(target)

        return self.config_dict


class Output:
    def __init__(self):
        self.targets = []
        self.config_paths = []
        self.config_dict = {}
        self.default_dict = {
            "minsize": 70,
            "maxsize": 200,
            "mpprimer": -3.5,
            "exception": [],
            "path": os.getcwd(),
            "intermediate": False,
            "qc_gene": ["rRNA"],
            "mfold": -3.0,
            "skip_download": False,
            "assemblylevel": ["all"],
            "skip_tree": False,
            "nolist": False,
            "offline": False,
            "ignore_qc": False,
            "mfethreshold": 90,
            "customdb": None,
            "blastseqs": 1000,
            "probe": False,
            "virus": False,
            "genbank": False,
            "evalue": 10,
            "nuc_identity": 0,
            "runmode": ["species"],
            "strains": []}

    def get_path(self):
        inpath = input(
                "Please specify a path to use as the working directory "
                "or hit return to use the current working directory:\n")
        if inpath:
            if os.path.isabs(inpath):
                path = inpath
            else:
                path = os.path.join(os.getcwd(), inpath)
        else:
            path = os.getcwd()
        return path

    def search_configfiles(self, path):
        print("Search in " + path)
        for root, dirs, files in os.walk(path):
            for file_name in files:
                if "config.json" == file_name:
                    file_path = os.path.join(root, file_name)
                    print("found:", file_path)
                    G.logger("> Found config file: " + file_path)
                    target = "/".join(root.split("/")[-2:-1])
                    abbr = H.abbrev(target)
                    primer_csv = os.path.join(
                        path, "Summary", target, abbr + "_primer.csv")
                    if not os.path.isfile(primer_csv):
                        self.targets.append(target)
                        self.config_paths.append(file_path)
                    else:
                        G.logger("> Skip run found results for " + target)

        return self.targets, self.config_paths

    def read_config(self, target, config_path):
        self.config_dict.update({target: {}})
        with open(config_path, "r") as f:
            for line in f:
                self.config_dict[target] = json.loads(line)

        # new 20.11.2018
        # make old config files compatible with the updated pipeline
        for key in self.default_dict:
            try:
                self.config_dict[target][key]
            except KeyError:
                self.config_dict[target].update({key: self.default_dict[key]})
                info1 = (
                    "Warning! new option " + str(key)
                    + " was not found in config file")
                info2 = (
                    "Used default value: " + str(key) + " = "
                    + str(self.default_dict[key]))
                print("\n" + info1 + "\n")
                print("\n" + info2 + "\n")
                G.logger(info1)
                G.logger(info2)

    def run_output(self):
        targets = input(
            "Search for config files for (a)ll or (s)elect targets:\n")
        path = self.get_path()
        if targets.lower() == "a" or targets.lower() == "all":
            self.search_configfiles(path)
        elif targets.lower() == "s" or targets.lower() == "select":
            targets = Input().targetinput()
            while (targets == "help" or targets == "" or targets is None):
                targets = Input().helpmessage(targets)

            self.targets = Input().parse_targets(targets)
            for target in self.targets:
                configpath = os.path.join(
                        path, target, "config", "config.json")
                if os.path.isfile(configpath):
                    self.config_paths.append(configpath)
                else:
                    info = (
                        "No configuration files found for "
                        + target + " in specified path")
                    print(info)
                    G.logger(info)
                    return self.run_output()
        else:
            print("no valid input: choose (a)ll or (s)elect")
            sys.exit()

        for index, target in enumerate(self.targets):
            config_path = self.config_paths[index]
            self.read_config(target, config_path)

        return self.config_dict

    def run_gui_output(self, targets, path):
        if targets is None:
            self.search_configfiles(path)
        else:
            self.targets = Input().parse_targets(targets)
            for target in self.targets:
                configpath = os.path.join(
                        path, target, "config", "config.json")
                if os.path.isfile(configpath):
                    self.config_paths.append(configpath)
                else:
                    self.config_paths.append('None')
                    info = (
                        "No configuration files found for "
                        + target + " in specified path")
                    print(info)
                    G.logger(info)

        for index, target in enumerate(self.targets):
            config_path = self.config_paths[index]
            if not config_path == 'None':
                self.read_config(target, config_path)
        return self.config_dict


if __name__ == "__main__":
    Input().initiate()
