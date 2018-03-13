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

outdir=OUT
refdir=OUTreference/BOOK
datadir=DATA

echo start TICLL-stuff

sort -t'#' -s -k4 -gr $refdir/rank.sorted > $outdir/chaintest.ranked
echo "start TICLL-chain"

$bindir/TICCL-chain -v $outdir/chaintest.ranked

if [ $? -ne 0 ]
then
    echo failed after TICLL-chain
    exit
fi

echo "checking chain results...."

sort $outdir/chaintest.ranked.chained > /tmp/chaintest.ranked.chained
diff /tmp/chaintest.ranked.chained $refdir/rank.chained >& /dev/null
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using diff /tmp/chaintest.ranked.chained $refdir/rank.chained"
    exit
fi

diff $outdir/chaintest.ranked.chained.debug $refdir/rank.chained.debug >& /dev/null
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using diff $outdir/chaintest.ranked.chained.debug $refdir/rank.chained.debug"
    exit
else
    echo OK!
fi

echo "and now caseless"

$bindir/TICCL-chain --caseless -o $outdir/caseless.ranked.chained $outdir/chaintest.ranked

if [ $? -ne 0 ]
then
    echo failed after TICLL-chain
    exit
fi

echo "checking caseless chain results...."

sort $outdir/caseless.ranked.chained > /tmp/caseless.ranked.chained
diff /tmp/caseless.ranked.chained $refdir/caseless.rank.chained >& /dev/null
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    echo "using diff /tmp/caseless.ranked.chained $refdir/caseless.rank.chained"
    exit
else
    echo OK!
fi
