FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
ENV BASH_ENV ~/.bashrc
WORKDIR /root

RUN apt update --yes && \
    apt install --yes --no-install-recommends \
        wget vim && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
RUN ./miniconda3/bin/conda init bash

RUN miniconda3/bin/conda install -n base -c conda-forge mamba
RUN miniconda3/bin/mamba create -c conda-forge -c bioconda -n snakemake snakemake
RUN miniconda3/bin/conda config --add channels conda-forge
RUN miniconda3/bin/conda config --add channels bioconda
RUN miniconda3/bin/conda config --add channels r
