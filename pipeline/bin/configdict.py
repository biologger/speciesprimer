#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 19:56:37 2019

@author: bio
"""
import os
import json
import csv

BASEPATH = "/home/pipeline"
DICTPATH = os.path.join(BASEPATH, "dictionaries")
abbrev_genus = os.path.join(DICTPATH, "genus_abbrev.csv")
primersett = os.path.join(BASEPATH, "p3parameters")
speclistfile = os.path.join(DICTPATH, "species_list.txt")

configdict = {
    "roary":{"default":True, "settings":{
        "p": None, "o": None, "e":None, "n":None, "i": None, "cd": None, 
        "g": None, "s": None, "t": None, "t": None, "ap": None, "z": None,
        "v": None, "y": None, "iv": None}}, 
    "specieslist":{"default":True, "list": [
        'Acidipropionibacterium acidipropionici', 'Acidipropionibacterium jensenii', 
        'Acidipropionibacterium thoenii', 'Acinetobacter baumannii', 
        'Acinetobacter bereziniae', 'Acinetobacter johnsonii', 
        'Acinetobacter parvus', 'Advenella kashmirensis', 
        'Aerococcus viridans', 'Agrococcus casei', 
        'Alcaligenes faecalis', 'Alkalibacterium kapii', 
        'Arthrobacter arilaitensis', 'Arthrobacter bergerei', 
        'Arthrobacter globiformis', 'Arthrobacter nicotianae', 
        'Bacillus altitudinis', 'Bacillus fordii', 
        'Bacillus galactosidilyticus', 'Bacillus pumilus', 
        'Bacillus subtilis', 'Bacteroides plebeius', 'Bavariicoccus seileri', 
        'Bifidobacterium adolescentis', 'Bifidobacterium animalis subsp. animalis',
        'Bifidobacterium animalis subsp. lactis', 'Bifidobacterium bifidum',
        'Bifidobacterium breve', 'Bifidobacterium longum', 
        'Bifidobacterium pseudolongum subsp. pseudolongum', 'Bifidobacterium thermophilum',
        'Brachybacterium alimentarium', 'Brachybacterium tyrofermentans', 
        'Brevibacillus parabrevis', 'Brevibacterium antiquum', 
        'Brevibacterium aurantiacum', 'Brevibacterium casei', 
        'Brevibacterium linens', 'Brevundimonas diminuta', 
        'Brochothrix thermosphacta', 'Carnobacterium divergens', 
        'Carnobacterium maltaromaticum', 'Carnobacterium mobile', 
        'Chryseobacterium bovis', 'Chryseobacterium ginsengisoli', 
        'Chryseobacterium haifense', 'Chryseobacterium joostei', 
        'Chryseobacterium shigense', 'Chryseobacterium vrystaatense', 
        'Citrobacter freundii', 'Citrobacter amalonaticus', 
        'Citrobacter farmeri', 'Citrobacter youngae', 'Citrobacter werkmanii',
        'Citrobacter portucalensis', 'Citrobacter pasteurii', 
        'Citrobacter braakii', 'Clostridium beijerinckii', 
        'Clostridium botulinum', 'Clostridium butyricum', 
        'Clostridium disporicum', 'Clostridium tyrobutyricum', 
        'Corynebacterium ammoniagenes', 'Corynebacterium bovis', 
        'Corynebacterium casei', 'Corynebacterium flavescens', 
        'Corynebacterium glutamicum', 'Corynebacterium propinquum', 
        'Corynebacterium stationis', 'Corynebacterium tuberculostearicum', 
        'Corynebacterium variabile', 'Cutibacterium acnes', 
        'Debaryomyces hansenii', 'Delftia acidovorans', 'Dietzia maris', 
        'Enterococcus casseliflavus', 'Enterococcus durans', 
        'Enterococcus faecalis', 'Enterococcus faecium', 
        'Enterococcus gallinarum', 'Enterococcus gilvus', 
        'Enterococcus hirae', 'Enterococcus italicus', 
        'Enterococcus malodoratus', 'Enterococcus pseudoavium', 
        'Escherichia coli', 'Exiguobacterium acetylicum', 
        'Exiguobacterium sibiricum', 'Facklamia tabacinasalis', 
        'Galactomyces candidum', 'Geotrichum candidum', 
        'Glutamicibacter arilaitensis', 'Glutamicibacter bergerei', 
        'Glutamicibacter nicotianae', 'Hafnia alvei', 'Halomonas venusta', 
        'Intestinibacter bartlettii', 'Janthinobacterium lividum', 
        'Jeotgalicoccus psychrophilus', 'Klebsiella aerogenes', 
        'Klebsiella oxytoca', 'Kluyvera intermedia', 'Kocuria atrinae', 
        'Kocuria carniphila', 'Kocuria kristinae', 'Kocuria palustris',
        'Kocuria rhizophila', 'Kocuria varians', 'Kurthia gibsonii',
        'Lactobacillus acidipiscis', 'Lactobacillus acidophilus',
        'Lactobacillus amylolyticus', 'Lactobacillus amylophilus',
        'Lactobacillus brevis', 'Lactobacillus buchneri',
        'Lactobacillus casei', 'Lactobacillus coryniformis',
        'Lactobacillus coryniformis subsp. coryniformis', 'Lactobacillus curvatus',
        'Lactobacillus delbrueckii', 'Lactobacillus delbrueckii subsp. bulgaricus',
        'Lactobacillus delbrueckii subsp. delbrueckii', 'Lactobacillus delbrueckii subsp. lactis',
        'Lactobacillus dextrinicus', 'Lactobacillus fermentum', 
        'Lactobacillus gallinarum', 'Lactobacillus gasseri', 
        'Lactobacillus graminis', 'Lactobacillus harbinensis',
        'Lactobacillus helveticus', 'Lactobacillus kefiri',
        'Lactobacillus nenjiangensis', 'Lactobacillus nodensis',
        'Lactobacillus oligofermentans', 'Lactobacillus parabrevis',
        'Lactobacillus parabuchneri', 'Lactobacillus paracasei',
        'Lactobacillus paracasei subsp. paracasei', 'Lactobacillus parafarraginis',
        'Lactobacillus parakefiri', 'Lactobacillus paraplantarum',
        'Lactobacillus pentosus', 'Lactobacillus perolens',
        'Lactobacillus plantarum', 'Lactobacillus plantarum subsp. argentoratensis',
        'Lactobacillus plantarum subsp. plantarum', 'Lactobacillus rhamnosus',
        'Lactobacillus sakei', 'Lactobacillus salivarius',
        'Lactobacillus satsumensis', 'Lactobacillus tucceti',
        'Lactobacillus vaccinostercus', 'Lactococcus garvieae',
        'Lactococcus lactis', 'Lactococcus lactis subsp. cremoris',
        'Lactococcus lactis subsp. lactis', 'Lactococcus raffinolactis',
        'Lelliottia amnigena', 'Leucobacter komagatae', 'Leuconostoc carnosum',
        'Leuconostoc citreum', 'Leuconostoc lactis', 'Leuconostoc mesenteroides',
        'Leuconostoc mesenteroides subsp. cremoris', 'Leuconostoc mesenteroides subsp. dextranicum',
        'Leuconostoc mesenteroides subsp. mesenteroides', 'Leuconostoc pseudomesenteroides',
        'Listeria innocua', 'Listeria monocytogenes', 'Luteococcus japonicus',
        'Macrococcus caseolyticus', 'Marinilactibacillus psychrotolerans',
        'Meiothermus silvanus', 'Microbacterium foliorum',
        'Microbacterium gubbeenense', 'Microbacterium lacticum',
        'Microbacterium oxydans', 'Micrococcus luteus',
        'Moraxella osloensis', 'Morganella morganii',
        'Morganella psychrotolerans', 'Mycobacterium vaccae',
        'Ornithinibacillus bavariensis', 'Paenibacillus lactis', 
        'Paludibacter propionicigenes', 'Pantoea agglomerans',
        'Pediococcus acidilactici', 'Pediococcus pentosaceus',
        'Pediococcus stilesii', 'Porphyromonas levii',
        'Propionibacterium acidipropionici', 'Propionibacterium acnes',
        'Propionibacterium freudenreichii', 'Propionibacterium freudenreichii subsp. freudenreichii',
        'Propionibacterium freudenreichii subsp. shermanii', 
        'Propionibacterium jensenii', 'Propionibacterium thoenii',
        'Proteus hauseri', 'Proteus vulgaris', 'Providencia alcalifaciens',
        'Providencia heimbachae', 'Providencia rettgeri',
        'Pseudoclavibacter helvolus', 'Pseudomonas asplenii', 
        'Pseudomonas fluorescens', 'Pseudomonas fragi',
        'Pseudomonas lundensis', 'Pseudomonas marginalis',
        'Pseudomonas putida', 'Psychrobacter aquimaris', 'Psychrobacter celer',
        'Psychrobacter faecalis', 'Psychrobacter immobilis', 
        'Psychrobacter namhaensis', 'Raoultella ornithinolytica', 
        'Raoultella planticola', 'Reyranella massiliensis', 
        'Rhodanobacter glycinis', 'Rothia mucilaginosa', 
        'Salmonella enterica subsp. enterica serovar Typhimurium', 
        'Serratia marcescens', 'Serratia proteamaculans', 
        'Serratia rubidaea', 'Sphingobacterium faecium', 
        'Staphylococcus aureus', 'Staphylococcus capitis', 
        'Staphylococcus caprae', 'Staphylococcus chromogenes', 
        'Staphylococcus cohnii', 'Staphylococcus epidermidis', 
        'Staphylococcus equorum', 'Staphylococcus equorum subsp. equorum', 
        'Staphylococcus equorum subsp. linens', 'Staphylococcus fleurettii', 
        'Staphylococcus haemolyticus', 'Staphylococcus hominis', 
        'Staphylococcus lentus', 'Staphylococcus sciuri', 
        'Staphylococcus sciuri subsp. sciuri', 'Staphylococcus succinus', 
        'Staphylococcus succinus subsp. casei', 'Staphylococcus vitulinus',
        'Staphylococcus warneri', 'Staphylococcus xylosus',
        'Stenotrophomonas indologenes', 'Stenotrophomonas maltophilia',
        'Stenotrophomonas rhizophila', 'Streptococcus infantarius subsp. infantarius',
        'Streptococcus lutetiensis', 'Streptococcus macedonicus', 
        'Streptococcus parauberis', 'Streptococcus pneumoniae', 
        'Streptococcus salivarius', 'Streptococcus thermophilus', 
        'Streptococcus uberis', 'Streptococcus vestibularis', 
        'Terrisporobacter mayombei', 'Tetragenococcus halophilus', 
        'Thermus thermophilus', 'Vagococcus fluvialis', 'Vagococcus lutrae',
        'Vibrio casei', 'Vibrio litoralis', 'Weissella cibaria', 
        'Weissella hellenica', 'Weissella paramesenteroides', 
        'Yarrowia lipolytica']},
    "primer3": {"default":True, "settings":{
            'PRIMER_PAIR_WT_PR_PENALTY': '1.0', 'PRIMER_MAX_SIZE': '26', 
            'PRIMER_MAX_NS_ACCEPTED': '0', 'PRIMER_PICK_ANYWAY': '0', 
            'PRIMER_OPT_GC_PERCENT': '50.0', 'PRIMER_QUALITY_RANGE_MIN': '0', 
            'PRIMER_INTERNAL_MAX_POLY_X': '5', 'PRIMER_OPT_TM': '60', 
            'PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE': '3', 'PRIMER_INTERNAL_WT_HAIRPIN_TH': '0.0',
            'PRIMER_TM_FORMULA': '1', 'PRIMER_PAIR_WT_COMPL_END': '0.0', 
            'PRIMER_INTERNAL_MIN_TM': '67.0', 'PRIMER_WT_TEMPLATE_MISPRIMING_TH': '0.0',
            'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': '20.00', 'PRIMER_INTERNAL_MAX_GC': '80.0',
            'PRIMER_PICK_RIGHT_PRIMER': '1', 'PRIMER_PAIR_MAX_COMPL_ANY': '8.00', 
            'PRIMER_WT_SELF_ANY': '0.0', 'PRIMER_WT_END_STABILITY': '0.0', 
            'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': '7', 'PRIMER_MAX_SELF_ANY': '8.00',
            'PRIMER_WT_POS_PENALTY': '0.0', 'PRIMER_INTERNAL_OPT_GC_PERCENT': '50.0',
            'PRIMER_PICK_LEFT_PRIMER': '1', 'PRIMER_INTERNAL_MAX_SELF_END_TH': '47.00',
            'PRIMER_TASK': 'generic', 'PRIMER_MASK_3P_DIRECTION': '0',
            'PRIMER_SALT_MONOVALENT': '50.0', 'PRIMER_INTERNAL_WT_TM_GT': '1.0',
            'PRIMER_WT_GC_PERCENT_GT': '0.0', 'PRIMER_MAX_GC': '70.0', 
            'PRIMER_INTERNAL_MIN_SIZE': '18', 'PRIMER_SALT_CORRECTIONS': '1', 
            'PRIMER_MASK_TEMPLATE': '1', 'PRIMER_INTERNAL_MAX_SELF_ANY': '12.00',
            'PRIMER_PAIR_WT_COMPL_ANY_TH': '0.0', 'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': '40.00',
            'PRIMER_MAX_END_GC': '3', 'PRIMER_INTERNAL_DNA_CONC': '50.0',
            'PRIMER_MAX_HAIRPIN_TH': '24.0', 'PRIMER_MAX_SELF_END_TH': '35.0', 
            'PRIMER_INTERNAL_MAX_NS_ACCEPTED': '0', 'PRIMER_WT_TM_LT': '1.0', 
            'PRIMER_INTERNAL_SALT_DIVALENT': '1.5', 'PRIMER_WT_LIBRARY_MISPRIMING': '0.0',
            'PRIMER_INTERNAL_WT_SELF_END_TH': '0.0', 'PRIMER_PAIR_WT_DIFF_TM': '0.0',
            'PRIMER_WT_SELF_END': '0.0', 'PRIMER_INTERNAL_WT_NUM_NS': '0.0', 
            'PRIMER_WT_SIZE_LT': '1.0', 'PRIMER_WT_GC_PERCENT_LT': '0.0', 
            'PRIMER_OUTSIDE_PENALTY': '0', 'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': '0', 
            'PRIMER_INTERNAL_MAX_LIBRARY_MISHYB': '12.00', 'PRIMER_MIN_SIZE': '18', 
            'PRIMER_INSIDE_PENALTY': '-1.0', 'PRIMER_SALT_DIVALENT': '1.5', 
            'PRIMER_INTERNAL_WT_SELF_END': '0.0', 'PRIMER_PAIR_MAX_COMPL_ANY_TH': '45.0',
            'PRIMER_INTERNAL_DNTP_CONC': '0.0', 'PRIMER_INTERNAL_WT_END_QUAL': '0.0', 
            'PRIMER_PAIR_MAX_DIFF_TM': '2', 'PRIMER_FIRST_BASE_INDEX': '1', 
            'PRIMER_INTERNAL_WT_SIZE_LT': '1.0', 'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': '0.0',
            'PRIMER_INTERNAL_WT_GC_PERCENT_LT': '0.0', 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': '1'
            , 'PRIMER_MASK_5P_DIRECTION': '1', 'PRIMER_MIN_GC': '30.0', 
            'PRIMER_MAX_POLY_X': '3', 'PRIMER_PAIR_MAX_COMPL_END': '3.00', 
            'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': '70.00', 'PRIMER_INTERNAL_MIN_QUALITY': '0',
            'PRIMER_DNTP_CONC': '0.6', 'PRIMER_LIBERAL_BASE': '1',
            'PRIMER_WT_END_QUAL': '0.0', 'PRIMER_INTERNAL_WT_GC_PERCENT_GT': '0.0', 
            'PRIMER_INTERNAL_MAX_SIZE': '27', 'PRIMER_MAX_TEMPLATE_MISPRIMING': '12.00',
            'PRIMER_GC_CLAMP': '0', 'PRIMER_MAX_TM': '62.0', 
            'PRIMER_INTERNAL_MAX_TM': '73.0', 'PRIMER_INTERNAL_SALT_MONOVALENT': '50.0',
            'PRIMER_MAX_END_STABILITY': '9.0', 'PRIMER_PAIR_WT_PRODUCT_TM_GT': '0.0', 
            'PRIMER_INTERNAL_OPT_SIZE': '20', 'PRIMER_PAIR_WT_PRODUCT_TM_LT': '0.0',
            'PRIMER_INTERNAL_WT_SEQ_QUAL': '0.0', 'PRIMER_WT_MASK_FAILURE_RATE': '0.0', 
            'PRIMER_DNA_CONC': '50.0', 'PRIMER_MAX_SELF_END': '3.00', 
            'PRIMER_INTERNAL_MAX_SELF_END': '12.00', 'PRIMER_MIN_END_QUALITY': '0',
            'PRIMER_LOWERCASE_MASKING': '0', 'PRIMER_MIN_QUALITY': '0', 
            'PRIMER_PAIR_MAX_COMPL_END_TH': '35.0', 'PRIMER_WT_TM_GT': '1.0',
            'PRIMER_INTERNAL_MAX_SELF_ANY_TH': '47.00', 'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH': '0.0',
            'PRIMER_INTERNAL_MAX_HAIRPIN_TH': '47.00', 'PRIMER_PAIR_WT_COMPL_ANY': '0.0', 
            'PRIMER_QUALITY_RANGE_MAX': '100', 'PRIMER_INTERNAL_WT_SELF_ANY_TH': '0.0',
            'PRIMER_OPT_SIZE': '20', 'PRIMER_INTERNAL_OPT_TM': '70.0', 
            'PRIMER_WT_SEQ_QUAL': '0.0', 'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': '24.00',
            'PRIMER_EXPLAIN_FLAG': '1', 'PRIMER_MAX_SELF_ANY_TH': '45.0', 
            'PRIMER_WT_SELF_ANY_TH': '0.0', 'PRIMER_PRODUCT_MIN_TM': '-1000000.0', 
            'PRIMER_WT_NUM_NS': '0.0', 'PRIMER_INTERNAL_WT_TM_LT': '1.0', 
            'PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE': '3', 'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': '0.0',
            'PRIMER_WT_HAIRPIN_TH': '1.0', 'PRIMER_PAIR_WT_COMPL_END_TH': '1.0', 
            'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING': '0.0', 'PRIMER_INTERNAL_WT_LIBRARY_MISHYB': '0.0',
            'PRIMER_WT_SELF_END_TH': '1.0', 'PRIMER_INTERNAL_MIN_GC': '20.0',
            'PRIMER_NUM_RETURN': '10', 'PRIMER_PRODUCT_MAX_TM': '90.0', 
            'PRIMER_MASK_FAILURE_RATE': '0.1', 'PRIMER_PAIR_WT_IO_PENALTY': '0.0',
            'PRIMER_WT_SIZE_GT': '1.0', 'PRIMER_MAX_LIBRARY_MISPRIMING': '12.00', 
            'PRIMER_PAIR_WT_LIBRARY_MISPRIMING': '0.0', 'PRIMER_WT_TEMPLATE_MISPRIMING': '0.0',
            'PRIMER_INTERNAL_WT_SIZE_GT': '1.0', 'PRIMER_MIN_TM': '57.0', 
            'PRIMER_INTERNAL_WT_SELF_ANY': '0.0', 'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': '4',
            'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': '0'
            }},
    "email": None,
    "genus_abbrev":{"default": True, "abbrev_dict":{
        'Pediococcus': 'Pd', 'Lactococcus': 'Lc', 'Propionibacterium': 'Pb', 
        'Escherichia': 'E', 'Leuconostoc': 'Ln', 'Clostridium': 'Cl', 
        'Hafnia': 'H', 'Streptococcus': 'Sc', 'Pseudomonas': 'Ps', 
        'Enterococcus': 'Ec'}}}






# allowed Roary options
"p" # INT number of threads [1]
"o" # STR clusters output filename [clustered_proteins]
"e" #     create a multiFASTA alignment of core genes using PRANK
"n" #     fast core gene alignment with MAFFT, use with -e
"i" #     minimum percentage identity for blastp [95]
"cd" # FLOAT percentage of isolates a gene must be in to be core [99]
"g" # INT maximum number of clusters [50000]
"s" #     dont split paralogs
"t" # INT translation table [11]
"ap" #    allow paralogs in core alignment
"z" #     dont delete intermediate files
"v" #     verbose output to STDOUT
"y" #     add gene inference information to spreadsheet, doesnt work with -e
"iv" # STR Change the MCL inflation value [1.5]

def get_primerdict():
    primerdict = {}
    with open(primersett) as f:
        for line in f:
            line = line.strip()
            if "=" in line:
                data = line.split("=")
                key = data[0]
                value = data[1]
                primerdict.update({key: value})
                
    print(primerdict)
    return primerdict

def get_specieslist():
    specieslist = []
    with open(speclistfile) as f:
        for line in f:
            line = line.strip()
            specieslist.append(line)
            
    print(specieslist)
    return specieslist

def get_genus_abbrev():
    abbrev = {}
    with open(abbrev_genus) as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            genus = row[0]
            abb = row[1]
            abbrev.update({genus: abb})
    
    print(abbrev)
    return abbrev
        
get_genus_abbrev()
