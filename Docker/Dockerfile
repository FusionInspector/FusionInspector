FROM ubuntu:20.04

MAINTAINER bhaas@broadinstitute.org

RUN  apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y gcc g++ perl python3 automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                       python3-setuptools python3-dev build-essential python3-distutils \
                                       unzip libbz2-dev  liblzma-dev gfortran libreadline-dev libcurl4-openssl-dev libx11-dev \
                       libxt-dev x11-common libcairo2-dev libpng-dev  libjpeg-dev pkg-config \
                       libxml2-dev  libssl-dev libcurl4-openssl-dev pbzip2 git \
                       libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev && \
    apt-get clean


## Python 3 stuff

RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN apt-get install -y python3-distutils

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py

RUN pip install requests igv-reports==1.8.0


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

ENV DATA /usr/local/data
RUN mkdir $DATA



## R installation:

WORKDIR $SRC

ENV R_VERSION=R-4.2.0

RUN curl https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && \
        tar xvf $R_VERSION.tar.gz && \
        cd $R_VERSION && \
        ./configure && make && make install


RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("cowplot")'
RUN R -e 'BiocManager::install("ranger")'


## perl lib installations

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install PerlIO::gzip
RUN cpanm install Set::IntervalTree
RUN cpanm install DB_File
RUN cpanm install URI::Escape
RUN cpanm install Carp::Assert
RUN cpanm install JSON::XS.pm


######################
## Tool installations:
######################


########
# Samtools

ENV SAMTOOLS_VERSION=1.9

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
   cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install


########
# Trinity

RUN apt-get install -y cmake apt-utils && apt-get clean

ENV TRINITY_VERSION=2.15.1

RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v${TRINITY_VERSION}/trinityrnaseq-v${TRINITY_VERSION}.FULL.tar.gz" && \
   cd $SRC && \
   wget $TRINITY_URL && \
   tar xvf trinityrnaseq-v${TRINITY_VERSION}.FULL.tar.gz && \
   cd trinityrnaseq-v${TRINITY_VERSION} && make


ENV TRINITY_HOME /usr/local/src/trinityrnaseq-v${TRINITY_VERSION}


## Bowtie2
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download -O bowtie2-2.3.3.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.3.1-linux-x86_64.zip && \
    mv bowtie2-2.3.3.1-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-2.3.3.1-linux-x86_64

## Jellyfish
WORKDIR $SRC
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz && \
    tar xvf jellyfish-2.2.7.tar.gz && \
    cd jellyfish-2.2.7/ && \
    ./configure && make && make install
            
## Salmon
WORKDIR $SRC
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz && \
    tar xvf Salmon-0.9.1_linux_x86_64.tar.gz && \
    ln -s $SRC/Salmon-latest_linux_x86_64/bin/salmon $BIN/.


## Picard tools

WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar
ENV PICARD_HOME $SRC



########
## Minimap2

WORKDIR $SRC
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - && \
    mv ./minimap2-2.26_x64-linux/minimap2 $BIN/

# K8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
    cp k8-0.2.4/k8-`uname -s` $BIN/k8



##############
## STAR

ENV STAR_VERSION=2.7.11a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


## Salmon
WORKDIR $SRC
ENV SALMON_VERSION=1.5.2
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    tar xvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -sf $SRC/salmon-${SALMON_VERSION}_linux_x86_64/bin/salmon $BIN/.



COPY sam_readname_cleaner.py $BIN/

RUN pip install numpy

# FusionInspector

ENV FI_VERSION=2.9.0
ENV FI_HASH=a43480df8dac6cfae0c01c2b636fd11de0d7bb98

RUN git clone --recursive https://github.com/FusionInspector/FusionInspector.git && \
    cd FusionInspector/ && \
    git checkout ${FI_HASH} && \
    git submodule init && git submodule update && \
    make && \
    mv * $BIN

