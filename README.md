# Classpoly

Hilber class polynomial computation in Z or modulo arbitrary integer P using CRT

Sources: 

- https://arxiv.org/pdf/0903.2785.pdf
- https://math.mit.edu/~drew/


## Installation

```bash
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

cd src/

tar -xzvf zn_poly-0.9.tar.gz
cd zn_poly-0.9
./configure
make && make install
cp include/* /usr/local/include/zn_poly/
cd ..

tar -xvf ff_poly_big_v1.2.7.tar
cd ff_poly_big_v1.2.7
make && make install
cp *.h /usr/local/include/ff_poly/
cd ..

tar -xvf classpoly_v1.0.2.tar
cd classpoly_v1.0.2
sed -e 's/^\(OBJECTS = \)/\1 prime.o /g' -i makefile
make
mkdir $HOME/temp
```

