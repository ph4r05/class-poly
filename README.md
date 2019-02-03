# Classpoly

Hilber class polynomial computation in Z or modulo arbitrary integer P using CRT

Sources: 

- https://arxiv.org/pdf/0903.2785.pdf
- https://math.mit.edu/~drew/

## Contents

This repository contains the following libraries:

- zn_poly-0.9  https://web.maths.unsw.edu.au/~davidharvey/code/zn_poly/releases/zn_poly-0.9.tar.gz
- ff_poly_big_v1.2.7  https://math.mit.edu/~drew/ff_poly_big_v1.2.7.tar
- classpoly_v1.0.2  https://math.mit.edu/~drew/classpoly_v1.0.2.tar

Original versions are in the branch `original`. Master contains slightly modified version.


## Modifications

The following modifications has been made:

- Logging to stderr (so output can be used for classpoly output)
- `CLASSPOLY_TEMP` defines the temporary directory to use. If undefined, goes by default to `$HOME/temp`
- `CLASSPOLY_PHI_FILES` defines directory with the PHI data files (downloaded from [SmallModPolys]). If undefined, default is `$HOME/phi_files`
- `CLASSPOLY_STDOUT=1` means the classpoly outputs the resulting polynomial to stdout instead of writing to the `$PWD/H_$D.txt`


## Docker build

```bash
$ docker build -t="classpoly" .
$ docker run -i -t classpoly
```

As described in the classpoly readme, you need to fetch polynomials from https://math.mit.edu/~drew/SmallModPolys
and extract it to $HOME/phi_files according to selected invariant (0=Hilbert class, `phi_j.tar`)

Inside the docker, you can fetch phi_files by calling `fetch-jinv.sh` for invariant=0.

If you have phi files already downloaded, mount the directory to the image and set env var correspondingly

```bash
$ docker build -t="classpoly" .
$ docker run -i --mount type=bind,src=/home/user/phi_files,dst=/usr/local/phi_files -t classpoly
$ export CLASSPOLY_PHI_FILES=/usr/local/phi_files
$ classpoly -19 0
Class polynomial for  D=-19 written to H_19.txt, degree 1, height 20 (20), size 0.000 MB (0.0s)

# For stdout:
$ CLASSPOLY_STDOUT=1 classpoly -19 0 2>/dev/null
I=0
D=-19
884736*X^0 +
1*X^1
```

## Build from sources

For build form the original sources:

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
cp ntutil.h /usr/local/include/
cd ..

tar -xvf classpoly_v1.0.2.tar
cd classpoly_v1.0.2
sed -e 's/^\(OBJECTS = \)/\1 prime.o /g' -i makefile
make
mkdir $HOME/temp
```


[SmallModPolys]: https://math.mit.edu/~drew/SmallModPolys