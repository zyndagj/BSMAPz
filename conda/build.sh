#!/bin/bash

samtools help 2>&1 | tail
make -j 4 bsmapz
make test
make DESTDIR=${PREFIX} install
