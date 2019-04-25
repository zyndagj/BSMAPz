#!/bin/bash

export CPATH=${PREFIX}/include
export LDFLAGS="-L${PREFIX}/lib"
export CPPFLAGS="-I$PREFIX/include"
export LIBRARY_PATH=${PREFIX}/lib

make -j2

make test

make DESTDIR=${PREFIX} install
