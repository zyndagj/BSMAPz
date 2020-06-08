#!/bin/bash

samtools help 2>&1 | tail

make bsmapz

make test

make DESTDIR=${PREFIX} install
