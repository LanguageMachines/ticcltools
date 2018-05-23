#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=book
fi

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

outdir=TESTRESULTS
refdir=OUTreference
testdir=TESTDATA
datadir=DATA

echo "start TICLL-rank"

$bindir/TICCL-rank -t 1 --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/ngram.ranked --debugfile $outdir/ngram.debug.ranked --artifrq 0 --clip 1 --skipcols=10,11 $refdir/ngram.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICLL-rank"
    exit
fi

echo "checking RANK results...."

diff $outdir/ngram.ranked $refdir/ngram.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using: diff $outdir/ngram.ranked $refdir/ngram.rank.sorted"
    exit
else
    echo OK!
fi
