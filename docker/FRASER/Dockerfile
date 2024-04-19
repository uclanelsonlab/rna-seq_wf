FROM --platform=linux/amd64 rocker/tidyverse:latest

# update and install required system libraries
RUN apt-get update -y && apt-get install -y \
    libglpk-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libbz2-dev

# install required R libraries
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install('FRASER')" && \
    R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')" && \
    R -e "BiocManager::install('org.Hs.eg.db')"