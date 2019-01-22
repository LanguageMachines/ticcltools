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

$bindir/TICCL-rank -t max --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/wv.ranked --debugfile $outdir/wv.debug.ranked --clip 5 --skipcols=10 $refdir/ngram.ldcalc --wordvec  /opensonar/SoNaRCurated/W2V/UctoNorm/SoNaR500.TRAIN.Dim100CBOW.vec --WVTEST

if [ $? -ne 0 ]
then
    echo "failed in TICLL-rank"
    exit
fi

echo "checking RANK clip5 results...."

diff $outdir/wv.ranked $refdir/wv.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using: diff $outdir/wv.ranked $refdir/wv.rank.sorted"
    exit
fi

sort $outdir/wv.debug.ranked > $outdir/wv.debug.rank.sorted
diff $outdir/wv.debug.rank.sorted $refdir/wv.debug.rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank debug results"
    echo "using: diff $outdir/wv.debug.rank.sorted $refdir/wv.debug.rank.sorted"
    exit
else
    echo OK!
fi
