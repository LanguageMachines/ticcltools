#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=unktest
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

outdir=OUT/$outsub
outdir=TESTRESULTS
refdir=TESTDATA
testdir=TESTDATA
datadir=DATA

echo start TICLL-unk

$bindir/TICCL-unk --corpus $datadir/INLandAspell.corpus --artifrq 100000000 --acro -o$outdir/unktest $testdir/unktest.tsv

if [ $? -ne 0 ]
then
    echo "failed in TICCL-unk"
    exit
fi
echo "checking UNK results...."
diff $outdir/unktest.punct $refdir/punct > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK punct results"
    echo "using: diff $outdir/unktest.punct $refdir/punct"
    exit
fi
diff $outdir/unktest.unk $refdir/unk > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK unk results"
    echo "using: diff $outdir/unktest.unk $refdir/unk"
    exit
fi
diff $outdir/unktest.clean $refdir/clean > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK clean results"
    echo "using: diff $outdir/unktest.clean $refdir/clean"
    exit
fi
diff $outdir/unktest.acro $refdir/acro > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK acro results"
    echo "using: diff $outdir/unktest.acro $refdir/acro"
    exit
fi

echo OK!
