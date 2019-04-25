#!/bin/bash

export CPATH=${PREFIX}/include
export LDFLAGS="-L${PREFIX}/lib"
export CPPFLAGS="-I$PREFIX/include"
export LIBRARY_PATH=${PREFIX}/lib

samtools help
make
make test
make DESTDIR=${PREFIX} install
