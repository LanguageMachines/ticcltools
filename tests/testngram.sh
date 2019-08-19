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

outdir=TESTRESULTS
refdir=OUTreference
testdir=TESTDATA
datadir=TRI

echo start TICLL-ldcalc

$bindir/TICCL-LDcalc --index $datadir/mre.indexNT --hash $datadir/mre.anahash --clean $datadir/mre.clean --LD 2 -t max --artifrq 100000000 -o $outdir/ngram.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICCL-LDcalc"
    exit
fi
echo "checking LDcalc results...."

LC_ALL=C sort $outdir/ngram.ldcalc  > /tmp/ldcalc

diff /tmp/ldcalc $refdir/ngram.ldcalc > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    echo "using: diff /tmp/ldcalc $refdir/ngram.ldcalc"
    exit
fi

LC_ALL=C sort $outdir/ngram.ldcalc.ambi  > /tmp/ldcalc.ambi
diff /tmp/ldcalc.ambi $refdir/ngram.ldcalc.ambi > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    echo "using: diff /tmp/ldcalc.ambi $refdir/ngram.ldcalc.ambi"
    exit
fi

LC_ALL=C sort $outdir/ngram.short.ldcalc  > /tmp/ldcalc.short
diff /tmp/ldcalc.short $refdir/ngram.ldcalc.short > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    echo "using: diff /tmp/ldcalc.short $refdir/ngram.ldcalc.short"
    exit
fi

echo "OK"
