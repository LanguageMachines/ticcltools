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

$bindir/TICCL-unk --background $datadir/INLandAspell.corpus --artifrq 100000000 --acro -o$outdir/smallunktest $testdir/smallunktest.tsv

if [ $? -ne 0 ]
then
    echo "failed in TICCL-unk"
    exit
fi
echo "checking UNK results...."
diff $outdir/smallunktest.punct $refdir/small.punct > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK punct results"
    echo "using: diff $outdir/smallunktest.punct $refdir/small.punct"
    exit
fi
diff $outdir/smallunktest.unk $refdir/small.unk > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK unk results"
    echo "using: diff $outdir/smallunktest.unk $refdir/small.unk"
    exit
fi
diff $outdir/smallunktest.clean $refdir/small.clean > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK clean results"
    echo "using: diff $outdir/smallunktest.clean $refdir/small.clean"
    exit
fi
diff $outdir/smallunktest.acro $refdir/small.acro > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK acro results"
    echo "using: diff $outdir/smallunktest.acro $refdir/small.acro"
    exit
fi

echo OK!
