#!/bin/bash

bindir=/home/sloot/usr/local/bin

if [ ! -d $bindir ]
then
   bindir=/exp/sloot/usr/local/bin
   if [ ! -d $bindir ]
   then
       echo "cannot find executables "
       exit
   fi
fi

outdir=OUT/TICCL
refdir=TESTDATA
datadir=DATA

echo "start TICLL-stats"

$bindir/TICCL-stats --ngram=1 $datadir/lines.txt -n -o $outdir/lines
$bindir/TICCL-stats --ngram=2 $datadir/lines.txt -n -o $outdir/lines
$bindir/TICCL-stats --ngram=3 $datadir/lines.txt -n -o $outdir/lines

if [ $? -ne 0 ]
then
    echo "failed in TICLL-stats"
    exit
fi

echo "checking TICCL-stats results...."
diff $outdir/lines.wordfreqlist.1.tsv $refdir/lines.1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-stats results"
    echo "using: diff $outdir/lines.wordfreqlist.1.tsv $refdir/lines.1"
    exit
fi
diff $outdir/lines.wordfreqlist.2.tsv $refdir/lines.2 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-stats results"
    echo "using: diff $outdir/lines.wordfreqlist.2.tsv $refdir/lines.2"
    exit
fi
diff $outdir/lines.wordfreqlist.3.tsv $refdir/lines.3 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-stats results"
    echo "using: diff $outdir/lines.wordfreqlist.3.tsv $refdir/lines.3"
    exit
fi

echo "done"
