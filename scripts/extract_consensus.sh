#!/bin/bash

awk 'function basename(file, a, n) {n = split(file, a, "/"); return a[n-1]}((FNR) % 2 == 1){print $0"_"basename(FILENAME)}((FNR) % 2 != 1){print $0}'
