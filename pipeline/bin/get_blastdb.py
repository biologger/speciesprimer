#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import hashlib
import wget
from datetime import datetime
import tarfile
import argparse

def commandline():
    parser = argparse.ArgumentParser(
        prog="Nucleotide collection (nt) database download",
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="Downloads the nucleotide collection (nt) database from NCBI",
        allow_abbrev=False)
    # define targetpath
    parser.add_argument(
        "-dbpath", "--dbpath", type=str, help="Path to save the DB")
    parser.add_argument(
        "-delete", "--delete", action="store_true",
        help="Delete the tar files after extraction to save Harddisk space")
    parser.add_argument(
        "-parts", "--parts", type=int, default=57,
        help="Number of parts of the nt database from ftp://ftp.ncbi.nlm.nih.gov/blast/db")

    args = parser.parse_args()
    return args

extractedendings = [
    ".nhd", ".nhi", ".nhr", ".nin", ".nnd",
    ".nni", ".nog", ".nsd", ".nsi", ".nsq"]

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
    with open(inputfile + ".md5", "r") as f:
        for line in f:
            md5sumcheck = line.split("  ")[0].strip()
            filenamecheck = line.split("  ")[1].strip()

    if inputfile == filenamecheck:
        print("\n\nmd5Checksum " + filenamecheck)
        filesum = md5Checksum(inputfile)
        if filesum == md5sumcheck:
            print("\nmd5checksum correct " + str(md5sumcheck) + "\n")
            return inputfile
        else:
            print("\nError MD5 checksum not correct")
            print("File: " + str(filesum))
            print("Expected: " + str(md5sumcheck))
            raise

def download_from_ftp(files, blastdb_dir, delete):
    os.chdir(blastdb_dir)
    baseurl = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db'
    for filename in files:
        check_extract = []
        if not os.path.isfile(filename):
            for end in extractedendings:
                if filename.split(".tar.gz")[0] + end in os.listdir("."):
                    check_extract.append(end)

            if check_extract == extractedendings:
                print("\nSkip download of " + filename + " found extracted files ")
            else:
                url = baseurl + "/" + filename
                print("\n" + str(datetime.now()))
                print("\nDownloading..." + filename)
                wget.download(url, filename)

                url = baseurl + "/" + filename + ".md5"
                print("\n" + str(datetime.now()))
                print("\nDownloading..." + filename + ".md5")
                wget.download(url, filename + ".md5")

                dbfile = check_md5(filename)
                extract_archives(dbfile, delete)

        else:
            dbfile = check_md5(filename)
            extract_archives(dbfile, delete)

def extract_archives(dbfile, delete):
    check_extract = []
    if dbfile is not None:
        for end in extractedendings:
            if not dbfile.split(".tar.gz")[0] + end in os.listdir("."):
                print("\nExtract archive " + dbfile)
                tar = tarfile.open(dbfile)
                tar.extractall()
                tar.close()
                check_extract.append(end)
            else:
                check_extract.append(end)
    if check_extract == extractedendings:
        print("Extracted " + dbfile)
        if delete == True:
            print("\nRemove archive " + dbfile)
            os.remove(dbfile)
            os.remove(dbfile + ".md5")


def generate_filelist(parts):
    files = []
    end = parts + 1
    for i in range(0, end):
        num = format(i, '02d')
        archive = "nt." + str(num) + ".tar.gz"
        files.append(archive)
    return files

def get_DB():
    print("Start Download of nucleotide collection (nt) database")
    args = commandline()
    if args.dbpath:
        if args.dbpath.endswith("/"):
            blastdb_dir = args.dbpath
        else:
            blastdb_dir = args.dbpath + "/"
    else:
        blastdb_dir = os.getcwd() + "/"

    filelist = generate_filelist(args.parts)
    download_from_ftp(filelist, blastdb_dir, args.delete)

    nt_nal_filepath = os.path.join(blastdb_dir, "nt.nal")

    if os.path.isfile(nt_nal_filepath):
        print("\nNucleotide collection database is ready")
    else:
        print("\nError nt.nal file for Nucleotide collection database is missing")

if __name__ == "__main__":
    get_DB()
