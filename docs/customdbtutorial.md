## custom DB tutorial

First we will download some sequences we want to use for our custom BLAST database with a handy tool from [kblin](https://github.com/kblin/ncbi-genome-download)

* Start a interactive terminal in your Docker container

		$ sudo docker exec -it {containername} bash

* In the container install ncbi-genome-download

		$ pip3 install ncbi-genome-download

* Download all complete genomes for the genus Lactobacillus in fasta format. For more details on how this works visit [https://github.com/kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)

		$ cd /blastdb
		$ ncbi-genome-download --genus Lactobacillus bacteria --assembly-level complete,chromosome --format fasta

After some time 414 Genomes were downloaded by ncbi-genome-download 

* The genomes can be found in a new refseq/bacteria directory

		$ cd refseq/bacteria

* Unzip all zipped fasta files

		$ gunzip -r *

* Copy the fasta files in the refseq directory

		$ cp **/*.fna /blastdb/refseq
		$ cd /blastdb/refseq

* Combine all files into one large fasta file (make sure no other files are in the refseq directory)

		$ cat *.fna > /blastdb/Lactobacillus_genomic.fas

* Create a BLAST database from the fasta file

		$ cd /blastdb
		$ makeblastdb -in Lactobacillus_genomic.fas -parse_seqids -title Lactobacillus_genomic -dbtype nucl -out lb_genomic

* This custom BLAST DB can now be used for primer design with SpeciesPrimer by specifying the path of the DB in the command line or in the GUI settings

		$ speciesprimer.py -customdb /blastdb/lb_genomic

* Maybe you want a more extensive database with more than just one genus you could add more sequences to your fasta file or combine different BLAST DB's with the blastdb_aliastool

		$ blastdb_aliastool -dblist "ref_prok_rep_genomes lb_genomic" -dbtype nucl -out repandlb -title "RefSeq Prok Rep and complete Lb genomes"
 
Now use in the command line -customdb /blastdb/repandlb.nal or with the GUI path to customBLASTDB: /blastdb/repandlb.nal and both databases are searched.
