#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import hashlib
import wget
import tarfile
import argparse
import logging
import time
import filecmp
import json
import shutil
import signal
import urllib.request
from html.parser import HTMLParser
from basicfunctions import GeneralFunctions as G

pipe_dir = os.path.dirname(os.path.abspath(__file__))
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')


class htmllinkparser(HTMLParser):
    StartTags = list()
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for name, value in attrs:
                if name == "href":
                    self.StartTags.append(value)


class config:
    baseend = [".nhr", ".nin", ".nnd", ".nni", ".nsq"]
    urldict = {
        "nt": {
            "base":"ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
            "http": "http://ftp.ncbi.nlm.nih.gov/blast/db/",
            "extend": [".nhd", ".nhi", ".nog"]},
        "ref_prok_rep_genomes": {
            "base": "ftp://ftp.ncbi.nlm.nih.gov/blast/db/",
            "http": "http://ftp.ncbi.nlm.nih.gov/blast/db/",
            "extend": [".nog", ".nsd", ".nsi"]},
        "test": {
            "base": "file:/blastdb/tmp/mockfiles/download/",
            "http": "file:/blastdb/tmp/mockfiles/download.html",
            "extend": [".nog", ".nsd", ".nsi"]}}

    def __init__(self, db, db_dir, delete, test):
        self.db = db
        self.db_dir = db_dir
        self.delete = delete
        if test:
            self.baseurl = self.urldict['test']['base']
            self.httpurl = self.urldict['test']['http']
            self.extract_end = self.baseend + self.urldict['test']['extend']
        else:
            self.baseurl = self.urldict[db]['base']
            self.httpurl = self.urldict[db]['http']
            self.extract_end = self.baseend + self.urldict[db]['extend']


def commandline():
    parser = argparse.ArgumentParser(
        prog="NCBI database download",
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="Downloads a pre-formatted database from NCBI",
        allow_abbrev=False)
    # define targetpath
    parser.add_argument(
        "-dbpath", "--dbpath", type=str, help="Path to save the DB")
    parser.add_argument(
        "-delete", "--delete", action="store_true", default=False,
        help="Delete the tar files after extraction to save Hard Disk space")
    parser.add_argument(
        "-db", "--database", type=str,
        help="Select nt or ref_prok_rep_genomes", default="nt",
        choices=["nt", "ref_prok_rep_genomes"])

    parser.add_argument('--test', action="store_true", help=argparse.SUPPRESS)

    return parser


def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
    return m.hexdigest()


def wget_download(filename, conf):
    archivename = filename.split(".md5")[0]
    url = conf.baseurl + "/" + archivename
    logger("> Downloading..." + archivename)
    wget.download(url, archivename)
    try:
        dbfile = check_md5(filename)
    except Exception:
        logger(
            "Problem with " + filename + ". Please run getblastdb"
            " again to check if files are missing")
        if os.path.isfile(filename):
            os.remove(filename)
        raise
    extract_archives(dbfile, conf)


def check_md5(inputfile):
    archivename = inputfile.split(".md5")[0]
    md5_file = inputfile
    with open(md5_file, "r") as f:
        for line in f:
            md5sumcheck = line.split("  ")[0].strip()
            filenamecheck = line.split("  ")[1].strip()
    if archivename == filenamecheck:
        logger("md5Checksum " + filenamecheck)
        filesum = md5Checksum(archivename)
        if filesum == md5sumcheck:
            logger("md5checksum correct " + str(md5sumcheck))
            return archivename

        logger("Error MD5 checksum not correct")
        logger("File: " + str(filesum))
        logger("Expected: " + str(md5sumcheck))
        raise Exception


def compare_md5_archive(inputfile, conf):
    archivename = inputfile.split(".md5")[0]
    with open(inputfile, "r") as f:
        for line in f:
            md5sumcheck = line.split("  ")[0].strip()
            filenamecheck = line.split("  ")[1].strip()
    if archivename == filenamecheck:
        logger("md5Checksum " + filenamecheck)
        filesum = md5Checksum(archivename)
        if filesum == md5sumcheck:
            logger("Skip download of " + archivename)
            dbfile = check_md5(inputfile)
            extract_archives(dbfile, conf)
        else:
            logger("Problem with " + archivename + " try to download it again")
            os.remove(archivename)
            try:
                wget_download(inputfile, conf)
            except Exception:
                logging.error("", exc_info=True)
            extract_archives(archivename, conf)


def compare_md5_files(filename, conf):
    old_file = os.path.join("md5_files", filename)
    if filecmp.cmp(filename, old_file) is True:
        logger("No change in " + filename + " skip download")
        os.remove(filename)
    else:
        os.remove(old_file)
        shutil.copy(filename, old_file)
        wget_download(filename, conf)

def get_extracted_endings(filename, conf):
    check_extract = []
    for end in conf.extract_end:
        if filename.split(".tar.gz.md5")[0] + end in os.listdir("."):
            check_extract.append(end)
    check_extract.sort()
    return check_extract


def handle_extracted_files(filename, file_path, conf):
    logger("Found extracted files, checking md5 file")
    if os.path.isfile(file_path):
        compare_md5_files(filename, conf)
    else:
        logger(" ".join([
                ">", "Skip download of", filename, "found extracted files"]))
        time.sleep(1)


def handle_old_md5file(filename, conf):
    archivename = filename.split(".md5")[0]
    if os.path.isfile(filename) and os.path.isfile(archivename):
        compare_md5_files(filename, conf)
        extract_archives(archivename, conf)
    else:
        wget_download(filename, conf)


def handle_md5_archive(filename, conf):
    archivename = filename.split(".md5")[0]
    if os.path.isfile(archivename):
        compare_md5_archive(filename, conf)
    else:
        wget_download(filename, conf)


def download_from_ftp(files, conf):
    os.chdir(conf.db_dir)
    for filename in files:
        file_path = os.path.join("md5_files", filename)
        try:
            check_extract = get_extracted_endings(filename, conf)
            conf.extract_end.sort()
            if check_extract == conf.extract_end:
                handle_extracted_files(
                    filename, file_path, conf)
            else:
                if os.path.isfile(file_path):
                    handle_old_md5file(
                            filename, conf)
                else:
                    handle_md5_archive(
                            filename, conf)

            if os.path.isfile(filename):
                os.rename(filename, file_path)

        except (KeyboardInterrupt, SystemExit):
            logging.error(
                "BLAST DB download was stopped while working on "
                + filename, exc_info=True)
            raise
        finally:
            for tmpfiles in os.listdir(conf.db_dir):
                if tmpfiles.endswith(".tmp"):
                    filepath = os.path.join(conf.db_dir, tmpfiles)
                    os.remove(filepath)


def extract_archives(dbfile, conf):
    check_extract = []
    extract_archive = []
    if dbfile is not None:
        for end in conf.extract_end:
            if not dbfile.split(".tar.gz")[0] + end in os.listdir("."):
                if dbfile not in extract_archive:
                    extract_archive.append(dbfile)
            else:
                check_extract.append(end)
        if len(extract_archive) > 0:
            for archive in extract_archive:
                logger("Extract archive " + dbfile)
                tar = tarfile.open(dbfile)
                tar.extractall()
                tar.close()
                check_extract = []
                for end in conf.extract_end:
                    if dbfile.split(".tar.gz")[0] + end in os.listdir("."):
                        check_extract.append(end)

    check_extract.sort()
    conf.extract_end.sort()
    if check_extract == conf.extract_end:
        logger("Extracted " + dbfile)
        if conf.delete is True:
            if os.path.isfile(dbfile):
                logger("Remove archive " + dbfile)
                os.remove(dbfile)


def get_md5files(conf):
    os.chdir(conf.db_dir)
    r = urllib.request.urlopen(conf.httpurl)
    parser = htmllinkparser()
    content = str(r.read())
    parser.feed(content)
    filelist = parser.StartTags
    for filename in filelist:
        if filename.startswith(conf.db + ".") and filename.endswith(".tar.gz.md5"):
            url = conf.baseurl + filename
            if not os.path.isfile(filename):
                logger("> Downloading..." + filename)
                wget.download(url, filename)


def get_filelist(conf):
    os.chdir(conf.db_dir)
    filelist = []
    for filename in os.listdir(conf.db_dir):
        if os.path.isfile(filename):
            if (
                filename.startswith(conf.db + ".")
                and filename.endswith(".tar.gz.md5")
            ):
                filelist.append(filename)
    filelist.sort()
    return filelist


def logger(string_to_log):
    print("\n" + string_to_log)
    logging.info(
        time.strftime(" %a, %d %b %Y %H:%M:%S: ", time.localtime())
        + string_to_log)


def exitatsigterm(signalNumber, frame):
    raise SystemExit('GUI stop')


def get_DB(mode=False):
    today = time.strftime("%Y_%m_%d", time.localtime())
    if mode == "auto":
        signal.signal(signal.SIGTERM, exitatsigterm)
        with open(tmp_db_path, 'r') as f:
            for line in f:
                tmp_db = json.loads(line)
        delete = tmp_db['BLAST_DB']['delete']
        db = tmp_db['BLAST_DB']['db']
        try:
            dbdir = tmp_db['BLAST_DB']['path']
            blastdb_dir = dbdir
        except KeyError:
            blastdb_dir = "/blastdb"
        try:
            if tmp_db['BLAST_DB']['test'] is True:
                test = True
        except KeyError:
            test = False

        logging.basicConfig(
            filename=os.path.join(
                "/", "primerdesign",
                "speciesprimer_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")

    else:
        parser = commandline()
        args = parser.parse_args()
        delete = args.delete
        test = args.test
        if args.dbpath:
            if args.dbpath.endswith("/"):
                blastdb_dir = args.dbpath
            else:
                blastdb_dir = args.dbpath + "/"
        else:
            blastdb_dir = os.getcwd() + "/"

        db = args.database

        logging.basicConfig(
            filename=os.path.join(
                blastdb_dir, "blastdb_download_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")

    G.create_directory(os.path.join(blastdb_dir, "md5_files"))

    conf = config(db, blastdb_dir, delete, test)

    logger("Start Download of NCBI " + db + " BLAST database")
    get_md5files(conf)
    filelist = get_filelist(conf)
    download_from_ftp(filelist, conf)

    nal_filepath = os.path.join(blastdb_dir, db + ".nal")
    if os.path.isfile(nal_filepath):
        logger("NCBI " + db + " BLAST database is ready")
    else:
        logger(
            "Error nal file for NCBI " + db + " BLAST database is missing")

if __name__ == "__main__":
    get_DB()
