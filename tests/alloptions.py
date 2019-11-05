#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:14:38 2019

@author: bio
"""

import itertools
# fulltest

customdb = [None]
probe = [True, False]
intermediate = [True, False]
skip_tree = [True, False]
blastdbv5 = [True, False]
exception = [None, "Lactobacillus_sakei"]
mfethreshold = [80, 90, 100]
mpprimer = [-3.5]
minsize = [70]
mfold = [-3.0]
target = ["Lactobacillus_curvatus"]
maxsize = [200]
qc_gene = [["rRNA"]]
nolist = [True, False]
assemblylevel = [["complete"]]
skip_download = [True, False]
ignore_qc = [True, False]
blastseqs = [100, 200, 500, 1000, 2000, 5000]
path = ["/home/primerdesign"]
offline = [True, False]



list_of_lists = [
        customdb, probe, intermediate, skip_tree, blastdbv5, exception,
        mfethreshold, mpprimer, minsize, mfold, target, maxsize,
        qc_gene, nolist, assemblylevel, skip_download, ignore_qc,
        blastseqs, path, offline]

r = len(list_of_lists)

combinations = list(itertools.combinations(list_of_lists, r))

print(combinations)
print(len(combinations))