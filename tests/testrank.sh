#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=book
fi

bindir=/home/sloot/usr/local/bin

outdir=TESTRESULTS
refdir=OUTreference
testdir=TESTDATA
datadir=DATA

echo "start TICLL-rank"

$bindir/TICCL-rank -t max --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/ngram.ranked --debugfile $outdir/ngram.debug.ranked --artifrq 100000000 --clip 1 --skipcols=10,11 $refdir/ngram.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICLL-rank"
    exit
fi

echo "checking RANK results...."

sort $outdir/ngram.ranked > /tmp/rank.sorted
diff /tmp/rank.sorted $refdir/ngram.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using: diff /tmp/rank.sorted $refdir/rank.sorted"
    exit
else
    echo OK!
fi
