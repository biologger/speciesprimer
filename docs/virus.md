## Virus example

### Limitations
  * Quality control for viral input genomes is __not supported__.
  * BLAST DB download for viruses is __not supported__

### Example
* Create a new directory for the BLAST DB

      $ mkdir /blastdb/virusdb/
      $ cd /blastdb/virusdb/

* Get a [pre-formatted BLAST DB](https://ftp.ncbi.nlm.nih.gov/blast/db/) or create a [custom BLAST DB](https://github.com/biologger/speciesprimer/tree/speciesprimerv2.2/docs/customdbtutorial.md)
      # download
      $ wget https://ftp.ncbi.nlm.nih.gov/blast/db/Betacoronavirus.tar.gz
      $ wget https://ftp.ncbi.nlm.nih.gov/blast/db/ref_viruses_rep_genomes.tar.gz
      # extract
      $ tar xvf ref_viruses_rep_genomes.tar.gz -C ref_viruses_rep_genomes
      $ tar xvf Betacoronavirus.tar.gz -C Betacoronavirus
      # combine
      $ blastdb_aliastool -dblist "Betacoronavirus/Betacoronavirus \
        ref_viruses_rep_genomes/ref_viruses_rep_genomes" -dbtype nucl \
        -out virusdb -title "Ref and betacoronavirus DB"

* Run the pipeline
      $ speciesprimer.py -t SARS-CoV-2 --virus --genbank --nuc_identity 90 \
      --customdb virusdb/virusdb --intermediate --nolist
