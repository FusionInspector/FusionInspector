FROM ubuntu:24.04

MAINTAINER bhaas@broadinstitute.org

RUN  apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y gcc g++ perl python3 automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlib1g zlib1g-dev  default-jre \
                       python3-setuptools python3-dev build-essential \
                                       unzip libbz2-dev  liblzma-dev gfortran libreadline-dev libcurl4-openssl-dev libx11-dev \
                       libxt-dev x11-common libcairo2-dev libpng-dev  libjpeg-dev pkg-config \
                       libxml2-dev  libssl-dev libcurl4-openssl-dev pbzip2 git \
                       libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev && \
    apt-get clean


## Python 3 stuff

RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py --break-system-packages


RUN pip install --break-system-packages requests igv-reports==1.8.0


## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

ENV DATA /usr/local/data
RUN mkdir $DATA



## R installation:

WORKDIR $SRC

ENV R_VERSION=R-4.4.2

RUN curl https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && \
        tar xvf $R_VERSION.tar.gz && \
        cd $R_VERSION && \
        ./configure && make && make install


RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("argparse")'

# Install additional system dependencies needed for tidyverse (arrow, ragg, feather)
# Install Arrow C++ libraries from Apache repository
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    cmake \
    libgit2-dev \
    libfontconfig1-dev \
    libwebp-dev \
    libutf8proc-dev \
    liblz4-dev \
    libzstd-dev \
    libsnappy-dev \
    libbrotli-dev \
    lsb-release \
    software-properties-common && \
    apt-get clean

# Install tidyverse
# Note: arrow is optional for tidyverse - core tidyverse packages don't require it
# tidyverse meta-package will install and load successfully even if arrow fails
ENV NOT_CRAN=true

# Install tidyverse (arrow may fail but that's OK - it's optional)
RUN R -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org", dependencies=TRUE)'

# Verify that library(tidyverse) works (it should load the core 8 packages)
RUN R -e 'library(tidyverse); cat("tidyverse version:", as.character(packageVersion("tidyverse")), "\n")' || \
    (echo "ERROR: library(tidyverse) failed to load" && exit 1)

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

ENV SAMTOOLS_VERSION=1.20

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
   cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install


########
# Trinity

RUN apt-get install -y cmake apt-utils autotools-dev && apt-get clean

ENV TRINITY_VERSION=2.15.2

RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v${TRINITY_VERSION}/trinityrnaseq-v${TRINITY_VERSION}.FULL.tar.gz" && \
   cd $SRC && \
   wget $TRINITY_URL && \
   tar xvf trinityrnaseq-v${TRINITY_VERSION}.FULL.tar.gz && \
   cd trinityrnaseq-v${TRINITY_VERSION} && \
   cp /usr/share/misc/config.guess /usr/share/misc/config.sub trinity-plugins/bamsifter/htslib/ && \
   sed -i 's/trinity_essentials: .*/trinity_essentials: jellyfish_target parafly_target/' Makefile && \
   make


ENV TRINITY_HOME /usr/local/src/trinityrnaseq-v${TRINITY_VERSION}


## Bowtie2
WORKDIR $SRC
ENV BOWTIE2_VERSION=2.5.4
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2* $BIN && \
    rm *.zip && \
    rm -r bowtie2-${BOWTIE2_VERSION}-linux-x86_64

## Jellyfish
WORKDIR $SRC
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz && \
    tar xvf jellyfish-2.2.7.tar.gz && \
    cd jellyfish-2.2.7/ && \
    ./configure && make && make install


## Picard tools

WORKDIR $SRC
ENV PICARD_VERSION=3.2.0
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar
ENV PICARD_HOME $SRC



########
## Minimap2

WORKDIR $SRC
#RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - && \
#    mv ./minimap2-2.26_x64-linux/minimap2 $BIN/

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - --no-same-owner && \
    mv ./minimap2-2.26_x64-linux/minimap2 $BIN/ && \
    mv ./minimap2-2.26_x64-linux/paftools.js $BIN/


# K8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && \
    cp k8-0.2.4/k8-`uname -s` $BIN/k8



##############
## STAR

ENV STAR_VERSION=2.7.11b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


## Salmon
WORKDIR $SRC
ENV SALMON_VERSION=1.10.3
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    tar xvf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
    ln -sf $SRC/salmon-${SALMON_VERSION}_linux_x86_64/bin/salmon $BIN/.



COPY sam_readname_cleaner.py $BIN/

RUN pip install --break-system-packages numpy

# FusionInspector

ENV FI_VERSION=2.11.0
ENV FI_HASH=836799ab22ce4289ff1d2431c934175201f0f6d6

RUN git clone --recursive https://github.com/FusionInspector/FusionInspector.git && \
    cd FusionInspector/ && \
    git checkout ${FI_HASH} && \
    git submodule init && git submodule update && \
    make && \
    mv * $BIN

