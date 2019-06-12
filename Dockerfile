FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
	apt-utils \
	nano \
	texlive-font-utils \
	wget \
	git \
	unzip \
	build-essential \
	gfortran \
	perl \
	bioperl \
	aragorn \
	prodigal \
	parallel \
	hmmer \
	infernal \
	barrnap \
	bedtools \
	cd-hit \
	mcl \
	gawk \
	cpanminus \
	prank \
	mafft \
	libdatetime-perl \
	libxml-simple-perl \
	libdigest-md5-perl \
	default-jre \
	emboss \
	python3-pip \
	python3-dev \
	python-pip


# create program directory
RUN mkdir /home/programs && mkdir /home/primerdesign
ENV PATH="/home/programs/:${PATH}"
ENV BLASTDB="/home/blastdb"

RUN cpanm -f Bio::Roary

# install python dependencies
COPY requirements.txt /
RUN pip3 install -r requirements.txt
RUN pip2 install psutil

# install prokka
RUN cd /home/programs && git clone https://github.com/tseemann/prokka.git \
&& prokka/bin/prokka --setupdb
ENV PATH="/home/programs/prokka/bin/:${PATH}"

# install blast-2.8.1+
RUN cd /home/programs && wget \
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.8.1+-x64-linux.tar.gz \
&& tar xvf ncbi-blast-2.8.1+-x64-linux.tar.gz
ENV PATH="/home/programs/ncbi-blast-2.8.1+/bin/:${PATH}"

# install primer3
RUN cd /home/programs && git clone https://github.com/primer3-org/primer3.git primer3 \
&& cd primer3/src && make && make test
ENV PATH="/home/programs/primer3/:${PATH}"

# libdg required by mfold
RUN cd /home/programs && wget -nv \
https://github.com/libgd/libgd/releases/download/gd-2.2.5/libgd-2.2.5.tar.gz \
&& tar xvf libgd-2.2.5.tar.gz && cd libgd-2.2.5 && ./configure && make && make install

# install mfold 3.6
RUN cd /home/programs && wget -nv \
http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz \
&& tar xvf mfold-3.6.tar.gz \
&& cd mfold-3.6 && ./configure && make && make install

# install MFEPrimer2.0
RUN cd /home/programs && wget -nv \
https://github.com/quwubin/MFEprimer/archive/v2.0.tar.gz \
&& tar xvf v2.0.tar.gz

ENV PATH="/home/programs/MFEprimer-2.0/:${PATH}"

#install tbl2asn
RUN cd /home/programs && wget -nv \
ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz \
&& gunzip linux.tbl2asn.gz

#install FastTreeMP
RUN cd /home/programs && wget -nv \
http://microbesonline.org/fasttree/FastTree.c \
&& gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o fasttree FastTree.c -lm

# Copy the directory contents into the docker directory
COPY pipeline /home/pipeline
ENV PATH="/home/pipeline/bin/:${PATH}"
ENV PATH="/home/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /home/pipeline/bin/*.py
RUN chmod +x /home/pipeline/ext-scripts/*.py

# Workdir
WORKDIR /home/primerdesign






