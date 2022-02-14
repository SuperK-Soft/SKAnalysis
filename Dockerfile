### Created by Dr. Benjamin Richards

### Download base image from repo
# starting image
FROM toolframework/core

### Run the following commands as super user (root):
USER root

# install 684MB software
# add openconnect to connect to VPN, needed to access SK svn repos
# we need the epel repo for openconnect, and we'll need sudo in a minute
RUN yum update -y \
&&  yum install -y \
    wget \
    tar \
    cmake \
    gcc-c++ \
    gcc \
    binutils \
    libX11-devel \
    libXpm-devel \
    libXft-devel \
    libXext-devel \
    libxml2-devel \
    libpng \
    libpng-devel \
    libjpeg-devel \
    graphviz-devel \
    mesa-libGL-devel \
    mesa-libGLU-devel \
    make \
    file \
    git \
    git-svn \
    bzip2-devel \
    cvs \
    automake \
    svn  \
    libtool \
    libxml2 \
    which \
    gsl-devel \
    gcc-gfortran \
    python36 \
    python36-libs \
    python36-devel \
    python36-pip \
    emacs \
    curl \
    curl-devel \
    ed \
    imake \
    makedepend \
    motif \
    motif-devel \
    krb5-libs \
    krb5-devel \
    libedit-devel \
    openldap \
    openldap-clients \
    openldap-devel \
    python-devel \
    openssl \
    openssl-devel \
    fftw \
    fftw-devel \
    mariadb-libs \
    mariadb-devel \
    postgresql-devel \
    libiodbc-devel \
    unixODBC-devel \
    openssl098e \
    nano \
    colordiff \
    man \
    gmp \
    gmp-devel \
    byacc \
    screen \
    gdb \
    compat-libf2c-34 \
    dialog \
    epel-release \
    openconnect \
    sudo \
    && yum clean all \
    && rm -rf /var/cache/yum

RUN groupadd sk && useradd --no-log-init --gid sk --uid 1000 --create-home sk_user
# make shell use bash by default
SHELL ["/bin/bash", "-c"]

# make a symlink for freetype required by CERN-LIB (details from cernlib website)
RUN ln -s /usr/include/freetype2/freetype /usr/include/freetype

# CONFIGURE VPN ### FIXME not sure this works....
###############
# to allow unprivilaged users to access vpn we need to do two things:
# first create a local tun network interface and allow sk_user to access it
# XXX if adding a user later to match that on the host we may wish to repeat this...
#RUN ip tuntap add vpn0 mode tun user sk_user  ##### h'mm, we can't do this yet
# as we need to map the /dev/net/tun from the host...

# the second thing is we still need root to configure the vpn....
# this involves having root run the /etc/vpnc/vpnc-script file, which gets invoked
# during connection. To do this without having to have the user we'll enable it to be run
# via sudo by users of the sk group without a password prompt.
RUN echo "%sk  ALL=(ALL) NOPASSWD:SETENV: /etc/vpnc/vpnc-script" >> /etc/sudoers

# STAGE TWO
###########
# Switch stages while we build. We'll copy over the final build products,
# dropping everything we don't need.
FROM base as build

# SETUP SKOFL SYMLINK
#####################
# for some reason skofl is linked to /usr/local on sukap machines,
# so mimic that in case anything looks for it there.
RUN mkdir -p /home/skofl/sklib_gcc4.8.5 \
&& ln -s /home/skofl/sklib_gcc4.8.5 /usr/local/sklib_gcc4.8.5

# INSTALL CERN-LIB
##################
RUN mkdir -p /home/skofl/sklib_gcc4.8.5/cern
# get a copy of cernlib 2005 (for Geant3)
# wget http://www-zeuthen.desy.de/linear_collider/cernlib/new/cernlib-2005-all-new.tgz
# the pre-build sourcefiles are only 37MB...
ADD hostfiles/cernlib-2005-all-new.tgz /home/skofl/sklib_gcc4.8.5/cern/

# remove the old patch file - actually we should be fine just overwriting it.
#RUN rm /home/skofl/sklib_gcc4.8.5/cern/cernlib.2005.corr.tgz
# get the new patch file
# wget http://www-zeuthen.desy.de/linear_collider/cernlib/new/cernlib.2005.corr.2014.04.17.tgz
# use COPY rather than ADD as it will be unpacked by the installer
#...few kB...
COPY hostfiles/cernlib.2005.corr.2014.04.17.tgz /home/skofl/sklib_gcc4.8.5/cern/cernlib.2005.corr.tgz

# same for the patch installation files
# wget http://www-zeuthen.desy.de/linear_collider/cernlib/new/cernlib.2005.install.2014.04.17.tgz
#... few more kB...
ADD hostfiles/cernlib.2005.install.2014.04.17.tgz /home/skofl/sklib_gcc4.8.5/cern/

# setup environment for building and build
WORKDIR /home/skofl/sklib_gcc4.8.5/cern
# somehow the installed dir is 700MB!
RUN . cernlib_env && ./Install_cernlib_and_lapack

# export CERN stuff for later installations
ENV CERN=/home/skofl/sklib_gcc4.8.5/cern
ENV CERN_LEVEL=2005
ENV CERN_ROOT=/home/skofl/sklib_gcc4.8.5/cern/2005
ENV PATH=$CERN_ROOT/bin:$PATH

# INSTALL ROOT
##############
RUN mkdir -p /home/skofl/sklib_gcc4.8.5/root_v5.28.00h/ && mkdir /home/skofl/sklib_gcc4.8.5/root_sauce/
#wget https://root.cern.ch/download/root_v5.28.00h.source.tar.gz
# initial unzipped sourcefiles are 120MB
ADD hostfiles/root_v5.28.00h.source.tar.gz /home/skofl/sklib_gcc4.8.5/root_sauce/

# ROOT 5 won't build out of the box on centos7, so we need to copied the modified files from sukap
# these were scp'd from the corresponding directories on sukap001
COPY hostfiles/rootfix/main.c /home/skofl/sklib_gcc4.8.5/root_sauce/root/build/rmkdepend/main.c
COPY hostfiles/rootfix/TUnixSystem.cxx /home/skofl/sklib_gcc4.8.5/root_sauce/root/core/unix/src/TUnixSystem.cxx
COPY hostfiles/rootfix/export.c /home/skofl/sklib_gcc4.8.5/root_sauce/root/graf2d/asimage/src/libAfterImage/export.c
COPY hostfiles/rootfix/import.c /home/skofl/sklib_gcc4.8.5/root_sauce/root/graf2d/asimage/src/libAfterImage/import.c
COPY hostfiles/rootfix/XrdOssAio.cc /home/skofl/sklib_gcc4.8.5/root_sauce/root/net/xrootd/src/xrootd/src/XrdOss/XrdOssAio.cc
COPY hostfiles/rootfix/GNUmakefile /home/skofl/sklib_gcc4.8.5/root_sauce/root/net/xrootd/src/xrootd/src/XrdCrypto/GNUmakefile

# ok, NOW we can compile and build
# configure options based on sukap001:/usr/local/sklib_gcc4.8.5/root_v5.28.00h/src/config.log
WORKDIR /home/skofl/sklib_gcc4.8.5/root_sauce/root/
# the built ROOT is 970MB
RUN ./configure --prefix=/usr/local/sklib_gcc4.8.5/root_v5.28.00h --libdir=/usr/local/sklib_gcc4.8.5/root_v5.28.00h/lib --incdir=/usr/local/sklib_gcc4.8.5/root_v5.28.00h/include --etcdir=/usr/local/sklib_gcc4.8.5/root_v5.28.00h/etc --enable-unuran --enable-roofit --enable-gdml --enable-minuit2 --with-cc=gcc --with-f77=gfortran --with-cxx=g++ \
&& make \
&& make install
# sukap retains the source files. Do we need them?
# this adds 810MB! We could at least remove a lot here by preventing the duplication,
# all it needs is putting the sourcefiles in the right location in the first place,
# it just means altering the zip file before inserting into the image
# (rename the top level directory from 'root' to 'src')
RUN mv /home/skofl/sklib_gcc4.8.5/root_sauce/root /home/skofl/sklib_gcc4.8.5/root_v5.28.00h/src \
&& rm -r /home/skofl/sklib_gcc4.8.5/root_sauce

# INSTALL NEUT
##############
#svn export https://kmcvs.icrr.u-tokyo.ac.jp/svn/rep/neut/tags/neut_5.4.0/
# i did this from sukap, then ran `tar -zcf neut_5.4.0.tgz neut_5.4.0` and scp'd the resulting zip locally
# but you can run it within the container (and omit the zip steps). Remember you need VPN started on the host.
# neut is 100MB pre-install
ADD hostfiles/neut_5.4.0.tgz /home/skofl/sklib_gcc4.8.5/
ENV NEUT_ROOT=/home/skofl/sklib_gcc4.8.5/neut_5.4.0
# make sure locale is set
# make sure fortran compiler is set
ENV LANG=C
ENV LC_ALL=C
ENV LC_MESSAGES=C
ENV Language=en_US.UTF-8
ENV FC=gfortran
WORKDIR $NEUT_ROOT/src/neutsmpl
# make sure thisroot.sh has been sourced, then build
# post-install neut is 360MB
RUN . /home/skofl/sklib_gcc4.8.5/root_v5.28.00h/bin/thisroot.sh \
&& /bin/csh Makeneutsmpl.csh

# this doesn't build zbsfns, which is required by SKOFL,
# but we can't do this until we have the SKOFL source files
# ($NEUT_ROOT/src/zbsfns/Makefile depends on $SKOFL_ROOT)
# As we don't have SKOFL yet, we'll need to delay this step until later,
# but note we will need to build it before we can compile the SKOFL libraries
#WORKDIR $NEUT_ROOT/src/zbsfns
#RUN . /home/skofl/sklib_gcc4.8.5/root_v5.28.00h/bin/thisroot.sh \
#&& /bin/csh make all \
#&& /bin/csh make install.library

# INSTALL GEANT4
#################
# first we need to build clhep, but before that we need a more up-to-date cmake
#wget https://github.com/Kitware/CMake/releases/download/v3.18.3/cmake-3.18.3-Linux-x86_64.tar.gz
# cmake pre-install is 113MB
ADD hostfiles/cmake-3.18.3-Linux-x86_64.tar.gz /usr/local/
# and of course moving it duplicates this 113MB, which could be removed by modifying the tar file
RUN mv /usr/local/cmake-3.18.3-Linux-x86_64 /usr/local/cmake-3.18.3
ENV PATH=/usr/local/cmake-3.18.3/bin:$PATH

# now we get can get CLHEP
#wget http://proj-clhep.web.cern.ch/proj-clhep/DISTRIBUTION/tarFiles/clhep-2.3.4.3.tgz
# the following extracts to 2.3.4.3, which will be the source directory
# clhep is only 5MB preinstall
ADD hostfiles/clhep-2.3.4.3.tgz /home/skofl/sklib_gcc4.8.5/
# make the install and build dirs
RUN mkdir /home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3 \
&&  mkdir /home/skofl/sklib_gcc4.8.5/clhep-build
WORKDIR /home/skofl/sklib_gcc4.8.5/clhep-build
# build, defaults std=c++11
# but somehow this installs 514MB!
RUN cmake -DCMAKE_INSTALL_PREFIX=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3 ../2.3.4.3/CLHEP \
&& cmake --build . --config RelWithDebInfo \
&& cmake --build . --target install
# clean up
RUN rm -r /home/skofl/sklib_gcc4.8.5/2.3.4.3/ \
&&  rm -r /home/skofl/sklib_gcc4.8.5/clhep-build
ENV CLHEP_BASE_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3
ENV CLHEP_INCLUDE_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3/include
ENV CLHEP_LIB_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3/lib

# get G4-10.03.p01
#wget https://www.jlab.org/12gev_phys/packages/sources/geant4/geant4.10.03.p01.tar.gz
# geant4 sourcefiles are 150MB
ADD hostfiles/geant4.10.03.p01.tar.gz /usr/local/sklib_gcc4.8.5/
# rename source dir as we want that name for the install dir, then make the build dir & install dir
# again moving adds another 150MB. We could avoid this by renaming the top directory in the tar file...
RUN mv /usr/local/sklib_gcc4.8.5/geant4.10.03.p01 /usr/local/sklib_gcc4.8.5/g4src \
&& mkdir /usr/local/sklib_gcc4.8.5/g4build && mkdir /usr/local/sklib_gcc4.8.5/geant4.10.03.p01
WORKDIR /usr/local/sklib_gcc4.8.5/g4build
# post-install geant4 is 2.25GB.
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local/sklib_gcc4.8.5/geant4.10.03.p01 -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_BUILD_CXXSTD=11 -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_ROOT_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3 ../g4src \
&& make -j2 \
&& make install
# cleanup
RUN rm -r /usr/local/sklib_gcc4.8.5/g4build \
&&  rm -r /usr/local/sklib_gcc4.8.5/g4src

# SKG4 on sukap actually uses a different geant4 installation entirely:
# the one in /usr/local/geant4.10/geant4.10.03.p02/
# notably this is patch 02 not 01, although it probably doesn't make a difference
# ***** DO NOT call G4ROOTsource.sh before building SKG4 ****
# ***** and we need to disable it in Make.sh too         ****
# invoke ./Make.sh only once at the start, then just call 'make'
# The only thing we need is add the missing LEND dataset (not part of p02) and point to it in bashrc
#ftp://gdo-nuclear.ucllnl.org/LEND_GND1.3/LEND_GND1.3_ENDF.BVII.1.tar.gz
# and the LEND files are 2.25GB just by themselves!
ADD hostfiles/LEND_GND1.3_ENDF.BVII.1.tar.gz /usr/local/sklib_gcc4.8.5/geant4.10.03.p01/share/data/

# random step for SKG4
#######################
RUN ln -s /usr/lib64/libXi.so.6 /usr/lib64/libXi.so

# random file for skdetsim-gd
#############################
RUN mkdir -p /home/atmpd/skdetsim/skdetsim-v13p90_for_atmnu_test/
# this random file is just a few bytes
COPY hostfiles/random.tbl.000 /home/atmpd/skdetsim/skdetsim-v13p90_for_atmnu_test/

#INSTALL PYTHON MODULES
########################
# these are used by the 2020 BDT for pure water neutron capture
USER sk_user
# we install 700MB of python stuff here
RUN pip3 install --user pybind11==2.6.2 \
                        xgboost==0.90 \
                        scikit-learn==0.23.2 \
                        numpy==1.19.5 \
                        scipy==1.5.4 \
                        joblib==0.14 \
                        matplotlib

# IIDA RELIC LIBRARIES
######################
# librelic libraries from iida-san's home directory
# these are often assumed to be present
USER root
RUN mkdir -p /home/iida/relic
# only a few hundred kB
ADD hostfiles/iida_relic.tar.gz /home/iida/

# SKOFL, ATMPD, skdetsim, SKG4 ... etc
######################################
# This will require interactively pulling using a user's SK credentials,
# so are not part of the dockerfile yet. All we'll do is copy in one of the setup scripts
# just a script, few kB
COPY hostfiles/bashenv_gcc4.8.5_skofl_19b+atmpd_19b /home/skofl/sklib_gcc4.8.5/bashenv_gcc4.8.5_skofl_19b+atmpd_19b

#
# For now these are also private repositories, but perhaps these will be public in the future...
#
## ToolDAQ
##########
## install the ToolDAQ framework
#WORKDIR /home/sk_user
#RUN git clone https://github.com/SuperK-UK-Soft/ToolAnalysis.git \
#&& cd ToolAnalysis && ./GetToolDAQ.sh --no_final
#
## SKROOT dictionaries
######################
## helper dictionaries to allow ROOT to recognise SKROOT classes
#RUN git clone https://github.com/SuperK-UK-Soft/skrootlibs.git \
#&& cd skrootlibs && make
#
## ROOT stl dictionaries
########################
## helper dictionaries to allow ROOT 5 to recognise more stl containers
#RUN git clone https://github.com:SuperK-UK-Soft/stllibs.git \
#&& cd stllibs && make

# CONVENIENCES
##############
# few kB
COPY hostfiles/bashrc /home/sk_user/.bashrc
RUN echo '. ~/.bashrc' >> /home/sk_user/.bash_profile
COPY hostfiles/screenrc /home/sk_user/.screenrc

# USER
######
# So that everything isn't owned by root, let's change the ownerships
# well of course this duplicates everything, so adds 6.68GB.
# we can get around this by adding the files with the right ownership to begin with
# see Dockerfile_v2
RUN chown -R sk_user:sk /home/skofl \
&&  chown -R sk_user:sk /home/iida \
&&  chown -R sk_user:sk /home/sk_user

# STAGE THREE
#############
# drop all the unnecessary duplications
# copy over everything we've installed
FROM base as install
COPY --from=build /home/skofl /home/skofl
COPY --from=build /home/iida /home/iida
COPY --from=build /home/sk_user /home/sk_user
COPY --from=build /usr/local/cmake-3.18.3 /usr/local/cmake-3.18.3
RUN ln -s /home/skofl/sklib_gcc4.8.5 /usr/local/sklib_gcc4.8.5

#  also need to copy over environmental variables
# CERNLIB
ENV CERN=/home/skofl/sklib_gcc4.8.5/cern
ENV CERN_LEVEL=2005
ENV CERN_ROOT=/home/skofl/sklib_gcc4.8.5/cern/2005
ENV PATH=$CERN_ROOT/bin:$PATH
# NEUT
ENV NEUT_ROOT=/home/skofl/sklib_gcc4.8.5/neut_5.4.0
ENV LANG=C
ENV LC_ALL=C
ENV LC_MESSAGES=C
ENV Language=en_US.UTF-8
ENV FC=gfortran
# GEANT4
ENV PATH=/usr/local/cmake-3.18.3/bin:$PATH
# CLHEP
ENV CLHEP_BASE_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3
ENV CLHEP_INCLUDE_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3/include
ENV CLHEP_LIB_DIR=/home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3/lib

USER sk_user
WORKDIR /home/sk_user

### Open terminal
CMD ["/bin/bash"]
