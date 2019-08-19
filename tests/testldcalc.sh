#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=ldcalctest
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
refdir=OUTreference/testldcalc
testdir=TESTDATA
datadir=DATA

echo start TICLL-ldcalc

$bindir/TICCL-LDcalc --alph=$datadir/nld.aspell.dict.clip20.lc.chars --index $refdir/id.indexNT --hash $refdir/anahash --clean $refdir/clean --LD 2 -t max --artifrq 100000000 -o $outdir/my.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICCL-LDcalc"
    exit
fi
echo "checking LDcalc results...."

LC_ALL=C sort $outdir/my.ldcalc  > /tmp/ldcalc

diff /tmp/ldcalc $refdir/ldcalc > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    echo "using: diff /tmp/ldcalc $refdir/ldcalc"
    exit
fi

echo "OK"
