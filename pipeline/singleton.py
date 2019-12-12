#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from basicfunctions import GeneralFunctions as G
from speciesprimer import CoreGenes

class Singletons(CoreGenes):
    
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.fasta_dir = os.path.join(self.results_dir, "fasta")
        self.all_core_path = os.path.join(self.pangenome_dir, "allcoregenes")
        self.multi_path = os.path.join(self.pangenome_dir, "multiannotated")
        self.single_dir = os.path.join(self.results_dir, "singletons")
        self.singleton = os.path.join(
                self.single_dir, "singleton_genes.csv")
        self.ffn_seqs = os.path.join(self.pangenome_dir, "ffn_sequences.csv")

    def get_singleton_genes(self):
        G.create_directory(self.single_dir)
        singleton_count = []
        filepath = os.path.join(
                self.pangenome_dir, "gene_presence_absence.csv")
        newtabledata = []
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            accessions = header[14:]
            genomes = len(accessions)
            for row in reader:
                data_row = []
                gene_name = row[0]
                number_isolates = int(row[3])
                number_sequences = int(row[4])
                loci = row[14:]
                if number_isolates == 1:
                    singleton_count.append(gene_name)
                    data_row.append(gene_name)
                    for locus in loci:              
                        if not locus == "":
                            print(locus)
                            data_row.append(locus)
                    newtabledata.append(data_row)


        G.csv_writer(self.singleton, newtabledata)
        return len(singleton_count)

    def get_fasta(self, locustags):

        def check_genename(gene):
            if "/" in gene:
                gene_name = "-".join(gene.split("/"))
            elif " " in gene:
                gene_name = "-".join(gene.split(" "))
            elif "'" in gene:
                gene_name = str(gene)[0:-1]
                gene = gene + "'"
            else:
                gene_name = gene

            return gene_name

        with open(self.singleton, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                gene = check_genename(row[0])
                for item in row[1:]:
                    name = locustags[item]["name"]
                    seq = locustags[item]["seq"]
                    gen_dir = os.path.join(self.single_dir, name)
                    G.create_directory(gen_dir)
                    outfile = os.path.join(gen_dir, gene + ".fasta")                    
                    with open(outfile, "w") as r:
                        record = SeqRecord(
                            Seq(seq),
                            name=item,
                            id='{}|{}|{}'.format(name, item, gene),
                            description="")
                        SeqIO.write(record, r, "fasta")

    def coregene_extract(self):
        info = "Run: core_gene_extract(" + self.target + ")"
        print(info)
        G.logger(info)
        self.get_singleton_genes()
        locustags = self.get_sequences_from_ffn()
        self.get_fasta(locustags)


def main():
    pass
    
    
if __name__ == "__main__":
    main()
