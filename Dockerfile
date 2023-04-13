FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]
ENV BASH_ENV ~/.bashrc
WORKDIR /src

RUN apt-get update -y && apt-get install -y bsdmainutils
RUN conda install -n base -c conda-forge mamba
RUN mamba create -c conda-forge -c bioconda -n snakemake snakemake
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN echo "source activate snakemake" >> ~/.bashrc
