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

# NO_BLAST.gi
BLAST databases can sometimes contain sequences with wrong taxonomic labels. This can prevent primer design by speciesprimer because no species specific sequences are identified or no specific primers can be found. Therefore such sequences can be excluded from the specificity evaluation.

Add the NCBI gene identifier(s) of a sequence to the NO_BLAST.gi file (one GI per line).
During the start of a new pipeline run the NO_BLAST.gi is read and copied to the config directory in the run. The GI's in this list are not used to evaluate the specificity of the template sequences or the primers.


