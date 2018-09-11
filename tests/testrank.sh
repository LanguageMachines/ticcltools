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

echo "start TICLL-rank clip=5"

$bindir/TICCL-rank -t max --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/ngram.c5.ranked --debugfile $outdir/ngram.debug5.ranked --artifrq 0 --clip 5 --skipcols=10,11 $refdir/ngram.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICLL-rank"
    exit
fi

echo "checking RANK clip5 results...."

diff $outdir/ngram.c5.ranked $refdir/ngram.c5.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using: diff $outdir/ngram.c5.ranked $refdir/ngram.c5.rank.sorted"
    exit
fi

sort $outdir/ngram.debug5.ranked > $outdir/ngram.debug5.rank.sorted
diff $outdir/ngram.debug5.rank.sorted $refdir/ngram.debug5.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank debug results"
    echo "using: diff $outdir/ngram.debug5.rank.sorted $refdir/ngram.debug5.rank.sorted"
    exit
else
    echo OK!
fi


echo "start TICLL-rank clip=1"

$bindir/TICCL-rank -t max --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/ngram.ranked --debugfile $outdir/ngram.debug.ranked --artifrq 0 --clip 1 --skipcols=10,11 $refdir/ngram.ldcalc

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
fi

sort $outdir/ngram.debug.ranked > $outdir/ngram.debug.ranked.sorted
diff $outdir/ngram.debug.ranked.sorted $refdir/ngram.debug.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank debug results"
    echo "using: diff $outdir/ngram.debug.ranked.sorted $refdir/ngram.debug.rank.sorted"
    exit
else
    echo OK!
fi
