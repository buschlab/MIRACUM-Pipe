#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit
  pwd -P
)

function install_texlive()
{
  apt-get install -y --no-install-recommends texlive texlive-lang-german texlive-latex-extra
  tlmgr option repository ftp://tug.org/historic/systems/texlive/2018/tlnet-final
  tlmgr init-usertree
  tlmgr install \
    breakurl \
    multirow \
    wrapfig \
    varwidth \
    threeparttable \
    preprint
}

function install_R()
{
  echo "deb http://cloud.r-project.org/bin/linux/debian bookworm-cran40/" >> /etc/apt/sources.list && \
  apt-key adv --keyserver keyserver.ubuntu.com --recv-key B8F25A8A73EACF41
  apt-get update && apt-get install -y --no-install-recommends -t bookworm-cran40 r-base-dev
  R CMD javareconf
}

# update packages
apt-get update & apt-get upgrade -y

# packages that are required for installation
# libcurl4-openssl-dev
apt-get install -y --no-install-recommends build-essential gcc-multilib libc-dev git-core cmake patch cmake ca-certificates \
  autoconf wget zip unzip zlib1g-dev libbz2-dev liblzma-dev libssl-dev libmariadb-dev tabix libmpfr-dev \
  libncurses5-dev libxml2-dev libcairo2-dev libxt-dev libgit2-dev libcurl4-gnutls-dev libncursesw5-dev libhts-dev libncurses5-dev \
  libharfbuzz-dev libfribidi-dev libtiff-dev\
  gfortran \
  default-jre \
  ant \
  perl-base \
  python3 python3-pysam python3-pip python3-numpy python3-scipy python3-matplotlib python3-reportlab python3-pandas python3-biopython python3-pyfaidx python3-pyvcf python3-setuptools python3-dev libpython3-all-dev python3-future python3-wheel \
  libsnappy-java \
  gnupg2 && \
  install_R && \
  install_texlive

# deps for VEP
  apt-get install -y --no-install-recommends libarchive-zip-perl libdbd-mysql-perl libdbi-perl libmodule-build-perl libwww-perl libbio-perl-perl

# python dependencies

  python3 -m pip install --upgrade pip
  python3 -m pip install -r /opt/MIRACUM-Pipe/debian/requirements.txt

# cleanup  

  apt-get -y autoremove
