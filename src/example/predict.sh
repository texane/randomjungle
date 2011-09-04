#!/usr/bin/env sh
./a.out \
--predict=o.jungle.xml \
--file=/tmp/fu.csv \
--delimiter=, \
--varnamesrow=0 \
--colselection=predict.cols \
--depvarname=Claim_Amount \
--outprefix=o \
--summary \
--verbose

