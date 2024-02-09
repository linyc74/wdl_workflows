FROM ghcr.io/atgenomix/runtime/base:1.5_22.04


RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    unzip \
    less \
    bwa=0.7.17-6 \
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
    cmake


# Link python3 to python
RUN ln -s /usr/bin/python3 /usr/bin/python \
 && python -m pip install pandas==1.3.5


RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt install --no-install-recommends r-base -y


ARG MINI_VERSION=4.11.0-0
ARG MINI_URL=https://github.com/conda-forge/miniforge/releases/download/${MINI_VERSION}/Mambaforge-${MINI_VERSION}-Linux-x86_64.sh

RUN curl --silent -L "${MINI_URL}" -o "mambaforge.sh" && \
    /bin/bash mambaforge.sh -b -p /opt/mambaforge/ && \
    rm mambaforge.sh

ENV PATH="/opt/mambaforge/bin:$PATH"

ENV PCGR_VERSION="1.1.0"
ENV PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock/"
ENV PLATFORM="linux"

RUN mamba create --name pcgr --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock && \
    mamba create --name pcgrr --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock && \
    mamba clean --all --force-pkgs-dirs --yes

ARG PCGR_ENV_NAME="pcgr"
ENV PATH="${PATH}:/opt/mambaforge/envs/${PCGR_ENV_NAME}/bin"
ENV CONDA_PREFIX="/opt/mambaforge/envs/${PCGR_ENV_NAME}"


RUN apt-get update \
 && apt-get install -y \
    unzip \
    tar \
    git \
    cpanminus \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libdbi-perl \
    libdbd-mysql-perl \
    libmysqlclient-dev

COPY netdb.h /usr/include/rpc/netdb.h

# use cpan for DBD:mysql because cpanm does not work
RUN cpanm DBI \
 && cpanm Archive::Zip \
 && cpan DBD::mysql \
 && cpanm Module::Build \
 && cpanm --force Bio::Root::Version \
 && cd /opt \
 && git clone -b release/110.1 --single-branch https://github.com/Ensembl/ensembl-vep.git \
 && cd ensembl-vep \
 && perl INSTALL.pl --AUTO ap --PLUGINS all --NO_UPDATE

ENV PATH="${PATH}:/opt/ensembl-vep"


RUN cd /usr/local/seqslab \
 && wget https://github.com/linyc74/omic/archive/refs/tags/v1.1.0.zip \
 && unzip v1.1.0.zip \
 && rm v1.1.0.zip \
 && mv omic-1.1.0 omic
