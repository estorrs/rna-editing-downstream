FROM python:3.6-jessie

# install packages
RUN apt-get update && apt-get install -y \
    vim

RUN pip install pytest sortedcontainers

# get conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p ./miniconda
ENV PATH="/usr/local/bin:/miniconda/bin:$PATH"

# add channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda install -y samtools

# set up working directory
COPY . /rna-editing-downstream
WORKDIR /rna-editing-downstream

CMD /bin/bash
