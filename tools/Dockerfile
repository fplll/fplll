## -*- docker-image-name: "fplll/fplll" -*-

FROM debian:buster
MAINTAINER Martin Albrecht <martin.albrecht@royalholloway.ac.uk>

ARG BRANCH=master

RUN apt update && \
    apt install -y build-essential libtool git autoconf libgmp-dev libmpfr-dev libqd-dev pkg-config && \
    apt clean && \
    git clone --branch $BRANCH https://github.com/fplll/fplll && \
    cd fplll && \
    autoreconf -i && \
    ./configure --prefix=/usr && \
    make install
    
