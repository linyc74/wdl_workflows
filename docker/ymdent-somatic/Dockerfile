FROM ghcr.io/atgenomix/runtime/base:1.5_22.04


RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    less \
    autoconf \
    automake \
    gcc \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    libc6-dev \
    libntirpc-dev \
    git \
    unzip \
    cmake \
    bwa=0.7.17-6 \
    bedtools=2.30.0+dfsg-2ubuntu0.1


# Remove old tabix, new tabix will be installed with htslib
RUN apt-get remove -y tabix


# Link python3 to python
RUN ln -s /usr/bin/python3 /usr/bin/python \
 && python -m pip install \
    pandas==1.3.5 \
    cutadapt==4.2


RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt install --no-install-recommends r-base -y


RUN cd /usr/local/seqslab \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && rm fastqc_v0.11.9.zip \
    && cd FastQC/ \
    && chmod 755 fastqc \
    && ln -s `pwd`/fastqc /usr/local/bin/fastqc


RUN cd /usr/local/seqslab \
    && curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz \
    && tar xvzf trim_galore.tar.gz \
    && rm trim_galore.tar.gz \
    && cd TrimGalore-0.6.6 \
    && ln -s `pwd`/trim_galore /usr/local/bin/trim_galore


RUN cd /usr/local/seqslab \
    && wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 \
    && wget https://github.com/samtools/samtools/releases/download/1.19.1/samtools-1.19.1.tar.bz2 \
    && wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
    && tar jxvf htslib-1.19.1.tar.bz2 \
    && tar jxvf samtools-1.19.1.tar.bz2 \
    && tar jxvf bcftools-1.19.tar.bz2 \
    && rm htslib-1.19.1.tar.bz2 samtools-1.19.1.tar.bz2 bcftools-1.19.tar.bz2 \
    && cd htslib-1.19.1 \
    && ./configure && make && make install && cd .. && rm -r htslib-1.19.1 \
    && cd bcftools-1.19 \
    && ./configure && make && make install && cd .. && rm -r bcftools-1.19 \
    && cd samtools-1.19.1 \
    && ./configure && make && make install && cd .. && rm -r samtools-1.19.1


RUN cd /usr/local/seqslab \
    && wget https://github.com/broadinstitute/gatk/releases/download/4.2.4.1/gatk-4.2.4.1.zip \
    && unzip gatk-4.2.4.1.zip \
    && rm gatk-4.2.4.1.zip \
    && cd gatk-4.2.4.1 \
    && ln -s `pwd`/gatk /usr/local/bin/gatk


RUN cd /usr/local/seqslab \
    && wget https://github.com/wwylab/MuSE/archive/refs/tags/v1.0-rc.tar.gz \
    && tar zxvf v1.0-rc.tar.gz \
    && rm v1.0-rc.tar.gz \
    && cd MuSE-1.0-rc \
    && make \
    && ln -s `pwd`/MuSE /usr/local/bin/muse


RUN cd /usr/local/seqslab \
    && wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz \
    && tar -xzf lofreq_star-2.1.5_linux-x86-64.tgz \
    && rm lofreq_star-2.1.5_linux-x86-64.tgz \
    && mv lofreq_star-2.1.5_linux-x86-64/bin/* /usr/local/bin


RUN cd /usr/local/seqslab \
    && wget https://github.com/linyc74/omic/archive/refs/tags/v1.2.0.zip \
    && unzip v1.2.0.zip \
    && rm v1.2.0.zip \
    && mv omic-1.2.0 omic
