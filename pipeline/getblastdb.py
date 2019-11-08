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

from basicfunctions import GeneralFunctions as G

BASEURLnt = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/'
BASEURLref = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/'

pipe_dir = os.path.dirname(os.path.abspath(__file__))
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')


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
        help="Delete the tar files after extraction to save Harddisk space")
    parser.add_argument(
        "-db", "--database", type=str,
        help="Select nt_v5 or ref_prok_rep_genomes", default="nt_v5",
        choices=["nt_v5", "ref_prok_rep_genomes"])
    args = parser.parse_args()
    return args


def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
    return m.hexdigest()


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
        else:
            logger("Error MD5 checksum not correct")
            logger("File: " + str(filesum))
            logger("Expected: " + str(md5sumcheck))
            raise


def compare_md5_archive(inputfile, delete, BASEURL, extractedendings):
    archivename = inputfile.split(".md5")[0]
    with open(inputfile, "r") as f:
        for line in f:
            md5sumcheck = line.split("  ")[0].strip()
            filenamecheck = line.split("  ")[1].strip()
    if archivename == filenamecheck:
        logger("md5Checksum " + filenamecheck)
        filesum = md5Checksum(archivename)
        if filesum == md5sumcheck:
            logger("No change in " + inputfile + " skip download")
            dbfile = check_md5(inputfile)
            extract_archives(dbfile, delete, extractedendings)
        else:
            os.rename(inputfile, os.path.join("md5_files", inputfile))
            url = BASEURL + "/" + archivename
            logger("Downloading..." + archivename)
            wget.download(url, archivename)
            dbfile = check_md5(inputfile)
            extract_archives(dbfile, delete, extractedendings)


def compare_md5_files(filename, delete, BASEURL, extractedendings):
    old_file = os.path.join("md5_files", filename)
    if filecmp.cmp(filename, old_file) is True:
        logger("No change in " + filename + " skip download")
        os.remove(filename)
    else:
        os.remove(old_file)
        os.rename(filename, os.path.join("md5_files", filename))
        archivename = filename.split(".md5")[0]
        url = BASEURL + "/" + archivename
        logger("Downloading..." + archivename)
        wget.download(url, archivename)
        dbfile = check_md5(filename)
        extract_archives(dbfile, delete, extractedendings)


def download_from_ftp(files, blastdb_dir, delete, BASEURL, extractedendings):
    os.chdir(blastdb_dir)
    for filename in files:
        check_extract = []
        try:
            for end in extractedendings:
                if filename.split(".tar.gz.md5")[0] + end in os.listdir("."):
                    check_extract.append(end)

            check_extract.sort()
            extractedendings.sort()
            if check_extract == extractedendings:
                logger("Found extracted files, checking md5 file")
                if os.path.isfile(os.path.join("md5_files", filename)):
                    compare_md5_files(filename, delete, BASEURL, extractedendings)
                else:
                    logger(
                        "> Skip download of " + filename
                        + " found extracted files ")
                    archivename = filename.split(".md5")[0]
                    time.sleep(1)
            else:
                file_path = os.path.join("md5_files", filename + ".md5")
                if os.path.isfile(file_path):
                    # maybe only the md5 file was downloaded
                    if os.path.isfile(filename):
                        compare_md5_files(filename, delete, BASEURL, extractedendings)
                    else:
                        url = BASEURL + "/" + filename
                        logger("> Downloading..." + filename)
                        wget.download(url, filename)
                        dbfile = check_md5(filename)
                        extract_archives(dbfile, delete, extractedendings)
                else:
                    archivename = filename.split(".md5")[0]
                    if os.path.isfile(archivename):
                        compare_md5_archive(filename, delete, BASEURL, extractedendings)
                    else:
                        url = BASEURL + "/" + archivename
                        logger("> Downloading..." + archivename)
                        wget.download(url, archivename)
                        dbfile = check_md5(filename)
                        extract_archives(dbfile, delete, extractedendings)
                        os.rename(
                            filename, os.path.join("md5_files", filename))
            if os.path.isfile(filename):
                os.rename(filename, os.path.join("md5_files", filename))
        except Exception as exc:
            msg = "> error while working on " + filename + " check logfile"
            print(msg)
            print(exc)
            logger(msg)
            logging.error("error while working on " + filename, exc_info=True)
            time.sleep(2)
            for files in os.listdir(blastdb_dir):
                if files.endswith(".tmp"):
                    filepath = os.path.join(blastdb_dir, files)
                    os.remove(filepath)

        except KeyboardInterrupt:
            logging.error(
                "KeyboardInterrupt while working on "
                + filename, exc_info=True)
            for files in os.listdir(blastdb_dir):
                if files.endswith(".tmp"):
                    filepath = os.path.join(blastdb_dir, files)
                    os.remove(filepath)
            raise
    for files in os.listdir(blastdb_dir):
        if files.endswith(".tmp"):
            filepath = os.path.join(blastdb_dir, files)
            os.remove(filepath)


def extract_archives(dbfile, delete, extractedendings):
    check_extract = []
    extract_archive = []
    if dbfile is not None:
        for end in extractedendings:
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
                for end in extractedendings:
                    if dbfile.split(".tar.gz")[0] + end in os.listdir("."):
                        check_extract.append(end)

    check_extract.sort()
    extractedendings.sort()
    if check_extract == extractedendings:
        logger("Extracted " + dbfile)
        if delete is True:
            logger("Remove archive " + dbfile)
            os.remove(dbfile)


def get_md5files(blastdb_dir, db, BASEURL):
    os.chdir(blastdb_dir)
    # -nc, --no-clobber skip downloads that would download to
    # existing files (overwriting them) can lead to problems
    # if single md5 files are manually deleted
    command = [
        "wget", "-nv", "-nc", "-r", "--no-parent", "--no-directories",
        "--tries", "4", "-A", db + ".*.tar.gz.md5",
        "-X", BASEURL+"FASTA," + BASEURL+"cloud," + BASEURL+"v5", BASEURL]
    G.run_subprocess(command)


def get_filelist(blastdb_dir, db):
    filelist = []
    for filename in os.listdir(blastdb_dir):
        if os.path.isfile(filename):
            if (
                filename.startswith(db + ".")
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


def get_DB(mode=False):
    today = time.strftime("%Y_%m_%d", time.localtime())
    if mode == "auto":
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

        logging.basicConfig(
            filename=os.path.join(
                "/", "primerdesign",
                "speciesprimer_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")

    else:
        args = commandline()
        delete = args.delete
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

    try:
        os.mkdir(os.path.join(blastdb_dir, "md5_files"))
    except Exception as exc:
        pass

    if db == "nt_v5":
        BASEURL = BASEURLnt
        extractedendings = [
            ".nhd", ".nhi", ".nhr", ".nin",
            ".nnd", ".nni", ".nog", ".nsq"]
    else:
        BASEURL = BASEURLref
        extractedendings = [
            ".nhr", ".nin", ".nnd", ".nni",
            ".nog", ".nsd", ".nsi", ".nsq"]

    logger("Start Download of NCBI " + db + " BLAST database")
    get_md5files(blastdb_dir, db, BASEURL)
    filelist = get_filelist(blastdb_dir, db)
    download_from_ftp(filelist, blastdb_dir, delete, BASEURL, extractedendings)

    nal_filepath = os.path.join(blastdb_dir, db + ".nal")
    if os.path.isfile(nal_filepath):
        logger("NCBI " + db + " BLAST database is ready")
    else:
        logger(
            "Error nal file for NCBI " + db + " BLAST database is missing")


if __name__ == "__main__":
    get_DB()
