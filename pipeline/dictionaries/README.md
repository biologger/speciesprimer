## genus_abbrev.csv

The genus abbreviations file is a .csv file with abbreviations for genus names.
The abbreviations make the names of the output files and primer pairs shorter.
e.g. Lb for Lactobacillus.

This comma separated list can be extended by adding new Genus names and Abbreviations.
e.g. 
Staphylococcus, Staph

This list is read during the initial stage of the pipeline runs and temporally stored as a dictionary.

## species_list.txt

The provided list contains bacterial species found in milk and cheese.
The species list is used to evaluate the specificity of the target sequences and to exclude unspecific primer pairs. It should contain the typical bacterial species you expect to be in the ecosystem you want to investigate with your PCR assay. Include also your target species, so the list can be reused for different targets, the pipeline filters this list at the beginning of a new run and writes the list of non-target species to the config file.

Create a species_list.txt file (one species name per line) or adapt it according to your needs. The pipeline reads the list during the initial steps of the run and these species names are then used to evaluate the BLAST results (NCBI pre-formatted nt database).

## 2.txids

The taxid list contains all NCBI taxonomic ids for bacteria (NCBI:txid2). From time to time NCBI may add new taxids if new species are added.

There is a script in the GUI version (Pipeline configuration), which allows to update this list if required.
This can also be done using python in the terminal of the docker container.

#### Example

$ python3
from Bio import Entrez
Entrez.email = "your_email_adress"
searchtaxid = Entrez.esearch(db="taxonomy", term="txid2[orgn]", retmax="500000")
taxidresult = Entrez.read(searchtaxid)
taxids = taxidresult["IdList"]
taxids.sort(key=lambda x: int(x))
filename = "2.txids"
with open(filename, "w") as f:
    for txid in taxids:
        f.write(txid + "\n")


As an alternative if the Entrez EDirect Utility and the BLAST+ command line application are installed on the host computer, the get_species_taxids.sh script from the BLAST+ program can be used to create this list.

__Example__

    $ get_species_taxids.sh -t 2 > $HOME/speciesprimer/pipeline/dictionaries/2.txids


