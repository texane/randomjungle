#!/usr/bin/env sh
./a.out \
--file=/home/texane/repo/kaggle/cpc/data/fu.csv \
--delimiter=, \
--ntree=100 \
--impmeasure=1 \
--varnamesrow=0 \
--colselection=train.cols \
--nrow=5000 \
--votes \
--verbose \
--outprefix=o \
--depvarname=Claim_Amount \
--summary \
--write=2
