# README NO_Blast
BLAST databases can sometimes contain sequences with wrong taxonomic labels. This can prevent primer design by speciesprimer because no species specific sequences are identified or no specific primers can be found. Therefore such sequences can be excluded from the specificity evaluation.

Add the NCBI gene identifier(s) of a sequence to the NO_BLAST.gi file (one GI per line).
During the start of a new pipeline run the NO_BLAST.gi is read and copied to the config directory in the run. The GI's in this list are not used to evaluate the specificity of the template sequences or the primers.
