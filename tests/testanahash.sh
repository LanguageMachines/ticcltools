#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=anahashtest
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
refdir=OUTreference/testanahash
testdir=TESTDATA
datadir=DATA

echo "start TICLL-anahash --list"

$bindir/TICCL-anahash --alph $datadir/nld.aspell.dict.clip20.lc.chars --list -o $outdir/analist $testdir/clean

if [ $? -ne 0 ]
then
    echo failed after TICCL-anahash
    exit
fi

echo "checking ANAHASH results...."
diff $outdir/analist.list $refdir/ok.list > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-anahash --list results"
    echo "using: diff $outdir/analist.list $refdir/ok.list"
    exit
else
    echo "OK"
fi
echo "start TICLL-anahash"

$bindir/TICCL-anahash --alph $datadir/nld.aspell.dict.clip20.lc.chars --artifrq 100000000 $testdir/clean

if [ $? -ne 0 ]
then
    echo failed after TICCL-anahash
    exit
fi

echo "checking ANAHASH results...."
LC_ALL=C sort $testdir/clean.corpusfoci > /tmp/foci
diff /tmp/foci $refdir/foci > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-anahash foci results"
    echo "using: diff /tmp/foci $refdir/foci"
    exit
fi

diff $testdir/clean.anahash $refdir/anahash > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-anahash foci results"
    echo "using: diff $testdir/clean.anahash $refdir/anahash"
    exit
else
    echo "OK"
fi
