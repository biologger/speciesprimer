#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import atexit
import logging
import subprocess
from mydaemon import Daemon


class Daemonize(Daemon):

    def run(self):
        RUNFILE = self.runfile
        self.command()
        # Delete the run file at exit. Maybe there will be no stop request.
        atexit.register(self.delrun)

        # Run while there is no stop request.
        n = 0
        while os.path.exists(RUNFILE):
            print('.', end='')
            n += 1

            if (n > 5):
                break

            time.sleep(1)


class GetBlastDB(Daemonize):

    pipe_dir = os.path.abspath(__file__).split("gui")[0]
    sys.path.append(pipe_dir)

    def command(self):
        import getblastdb
        getblastdb.get_DB(mode="auto")


class SpeciesPrimer(Daemonize):

    pipe_dir = os.path.abspath(__file__).split("gui")[0]
    sys.path.append(pipe_dir)

    def command(self):
        import speciesprimer
        speciesprimer.main(mode="auto")


class Update_prokDB(Daemonize):

    def command(self):
        today = time.strftime("%Y_%m_%d", time.localtime())
        logging.basicConfig(
            filename=os.path.join(
                "/", "primerdesign",
                "speciesprimer_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")
        blastdb_dir = "/blastdb"
        os.chdir(blastdb_dir)
        command = [
            "update_blastdb.pl", "--passive", "--decompress",
            "ref_prok_rep_genomes"]
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    print(output)
                    logging.info(
                        time.strftime("%d %b %Y %H:%M:%S: ", time.localtime())
                        + str(output).strip())
                else:
                    break
        while process.poll() is None:
            check_output()


class Update_ntDB(Daemonize):

    def command(self):
        today = time.strftime("%Y_%m_%d", time.localtime())
        logging.basicConfig(
            filename=os.path.join(
                "/", "primerdesign",
                "speciesprimer_" + today + ".log"),
            level=logging.DEBUG, format="%(message)s")
        blastdb_dir = "/blastdb"
        os.chdir(blastdb_dir)
        command = ["update_blastdb.pl", "--passive", "--decompress", "nt"]
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        def check_output():
            while True:
                output = process.stdout.readline().decode().strip()
                if output:
                    print(output)
                    logging.info(
                        time.strftime("%d %b %Y %H:%M:%S: ", time.localtime())
                        + str(output).strip())
                else:
                    break
        while process.poll() is None:
            check_output()
