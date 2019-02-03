# Multistage docker build, requires docker 17.05

# Base image is an argument - debian or ubuntu.
ARG BASE_IMAGE=debian:latest

# Builder stage
# https://docs.docker.com/develop/develop-images/dockerfile_best-practices/
# https://docs.docker.com/engine/reference/builder/
# https://medium.com/@tonistiigi/advanced-multi-stage-build-patterns-6f741b852fae
FROM ${BASE_IMAGE} AS base

RUN set -ex && \
    apt-get update && \
    apt-get --no-install-recommends --yes install \
        ca-certificates \
        cmake \
        g++ \
        make \
        pkg-config \
        git \
        curl \
        libtool-bin \
        autoconf \
        automake \
        bzip2 \
        xsltproc \
        gperf \
        unzip \
        rsync \
        python \
        libntl-dev \
        libgmp-dev

# Building class poly
FROM base AS builder
WORKDIR /usr/local

ENV PROJECT_DIR /usr/local/class-poly
ENV ZN_POLY zn_poly-0.9
ENV FF_POLY_BIG ff_poly_big_v1.2.7
ENV CLASSPOLY classpoly_v1.0.2

# Main GitHub repo
RUN set -ex \
    && git clone https://github.com/ph4r05/class-poly

# zn_poly
WORKDIR $PROJECT_DIR/src/$ZN_POLY
RUN set -ex \
    && ./configure \
    && make \
    && make install \
    && cp include/* /usr/local/include/zn_poly/

# ff_poly_big_v1.2.7
WORKDIR $PROJECT_DIR/src/$FF_POLY_BIG
RUN set -ex \
    && make \
    && make install \
    && cp *.h /usr/local/include/ff_poly/ \
    && cp ntutil.h /usr/local/include/

# classpoly_v1.0.2
WORKDIR $PROJECT_DIR/src/$CLASSPOLY
RUN set -ex \
    && sed -e 's/^\(OBJECTS = \)/\1 prime.o /g' -i makefile \
    && make \
    && mkdir $HOME/temp \
    && cp classpoly /usr/local/bin/


WORKDIR /usr/local
########################################################################################################################
## OpenSSL
#ARG OPENSSL_VERSION=1.1.0j
#ARG OPENSSL_HASH=31bec6c203ce1a8e93d5994f4ed304c63ccf07676118b6634edded12ad1b3246
#RUN set -ex \
#    && curl -s -O https://www.openssl.org/source/openssl-${OPENSSL_VERSION}.tar.gz \
#    && echo "${OPENSSL_HASH}  openssl-${OPENSSL_VERSION}.tar.gz" | sha256sum -c \
#    && tar -xzf openssl-${OPENSSL_VERSION}.tar.gz \
#    && cd openssl-${OPENSSL_VERSION} \
#    && ./Configure linux-x86_64 no-shared --static -fPIC \
#    && make build_generated \
#    && make libcrypto.a \
#    && make install
#ENV OPENSSL_ROOT_DIR=/usr/local/openssl-${OPENSSL_VERSION}
#
## ZMQ
#ARG ZMQ_VERSION=v4.2.5
#ARG ZMQ_HASH=d062edd8c142384792955796329baf1e5a3377cd
#RUN set -ex \
#    && git clone https://github.com/zeromq/libzmq.git -b ${ZMQ_VERSION} \
#    && cd libzmq \
#    && test `git rev-parse HEAD` = ${ZMQ_HASH} || exit 1 \
#    && ./autogen.sh \
#    && CFLAGS="-fPIC" CXXFLAGS="-fPIC" ./configure --enable-static --disable-shared \
#    && make \
#    && make install \
#    && ldconfig


#WORKDIR /src
#COPY . .

#ARG NPROC
#RUN set -ex && \
#    git submodule init && git submodule update && \
#    rm -rf build && \
#    if [ -z "$NPROC" ] ; \
#    then make -j$(nproc) release-static ; \
#    else make -j$NPROC release-static ; \
#    fi

## runtime stage
#FROM ubuntu:16.04
#
#RUN set -ex && \
#    apt-get update && \
#    apt-get --no-install-recommends --yes install ca-certificates && \
#    apt-get clean && \
#    rm -rf /var/lib/apt
#COPY --from=builder /src/build/release/bin /usr/local/bin/

# Contains data
#VOLUME /root/.test