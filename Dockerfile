FROM rocker/r-ver:4.3.2

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# gawk for taxon-tools
# gcc through libtool for treePL
# zlib1g-dev for R package XVector
# libxml2-dev for R package XML
# libudunits2-dev for R package units
# libgdal-dev for R package sf
# libzmq3-dev for R package rzmq -> clustermq
# libmagick++-dev for R package magick -> phytools
# python-dev-is-python3 for biopython -> superCRUNCH
# git-lfs for gittargets
# libarchive-dev for archive
# libharfbuzz-dev, libfribidi-dev for R package textshaping
# librdf0-dev for redland -> R package deposits
# cmake -> R package nanonext

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libnlopt-dev \
    libnlopt0 \
    libcolpack-dev \
    make \
    libomp-dev \
    build-essential \
    autoconf \
    autotools-dev \
    automake \
    libtool \
    zlib1g-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libzmq3-dev \
    libmagick++-dev \
    mafft \
    ncbi-blast+ \
    fastp \
    time \
    parallel \
    python-dev-is-python3 \
    curl \
    fasttree \
    gawk \
    cd-hit \
    git-lfs \
    libarchive-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    librdf0-dev \
    libgit2-dev \
    cmake \
    wget \
    libglpk40 \
    libarchive13 \
    libzmq5 \
    liblzma-dev \
    libbz2-dev \
    libsecret-1-0 \
  && apt-get clean

########################
### python libraries ###
########################

# biopython for superCRUNCH

RUN curl https://bootstrap.pypa.io/get-pip.py | python \
  && pip install biopython

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### IQ Tree v2 ###
ENV APP_NAME=iqtree
ENV IQT_VERSION=2.2.2.7
RUN wget https://github.com/iqtree/iqtree2/releases/download/v$IQT_VERSION/iqtree-$IQT_VERSION-Linux.tar.gz \
  && tar xf $APP_NAME-$IQT_VERSION-Linux.tar.gz \
  && rm $APP_NAME-$IQT_VERSION-Linux.tar.gz \
  && mv $APP_NAME-$IQT_VERSION-Linux/bin/iqtree2 /usr/local/bin/

### trimAL ###
ENV APP_NAME=trimal
ENV TRIMAL_VERSION=492b6d9455ec2d0b19a420bcfc4eea8adc88e3ea
RUN git clone https://github.com/scapella/$APP_NAME.git && \
	cd $APP_NAME && \
  git checkout $TRIMAL_VERSION && \
  cd source && \
	make && \
	cp trimal /usr/local/bin

### pandoc ###

ENV APP_NAME=pandoc
ENV PANDOC_VERSION=3.1.6.1
RUN wget https://github.com/jgm/pandoc/releases/download/$PANDOC_VERSION/pandoc-$PANDOC_VERSION-1-amd64.deb \
  && dpkg -i pandoc-$PANDOC_VERSION-1-amd64.deb \
  && rm pandoc-$PANDOC_VERSION-1-amd64.deb

### taxon-tools ###
# needs pandoc
ENV APP_NAME=taxon-tools
ENV TAXONTOOLS_VERSION=5bc8c1f7b51ed51d22773f6b16d75ccec3ae6dec
RUN git clone https://github.com/camwebb/$APP_NAME.git && \
	cd $APP_NAME && \
  git checkout $TAXONTOOLS_VERSION && \
	make check && \
	make install

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir /renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir /tmp/project

COPY ./renv.lock /tmp/project

WORKDIR /tmp/project

# Restore, but don't use cache
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

# Install packages used by VS code, but don't snapshot
RUN Rscript -e 'pkg_ignore <- c("jsonlite", "rlang", "languageserver", "reprex", "mdlincoln/docthis"); renv::install(pkg_ignore); renv::settings$ignored.packages(pkg_ignore, persist = TRUE)'

WORKDIR /home/