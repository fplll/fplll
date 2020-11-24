## -*- docker-image-name: "fplll/fplll" -*-

FROM debian:buster
MAINTAINER Martin Albrecht <fplll-devel@googlegroups.com>

ARG BRANCH=master
ARG JOBS=2
ARG CXXFLAGS="-O2 -march=x86-64"
ARG CFLAGS="-O2 -march=x86-64"

SHELL ["/bin/bash", "-c"]
ENTRYPOINT /bin/bash

RUN apt update && \
    apt install -y build-essential libtool git autoconf libgmp-dev libmpfr-dev libqd-dev pkg-config && \
    apt clean && \
    git clone --branch $BRANCH https://github.com/fplll/fplll && \
    cd fplll && \
    autoreconf -i && \
    CFLAGS=$CFLAGS CXXFLAGS=$CXXFLAGS ./configure --disable-static --prefix=/usr && \
    make -j $JOBS install && \
    cd .. && \
    rm -rf fplll 

    
