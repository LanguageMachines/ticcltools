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
testdata=TESTDATA


echo "start TICCL-chainclean"

$bindir/TICCL-chainclean --lexicon $testdata/chain.clean -o /tmp/chained.clean $testdata/morse.chained $1 $2 $3

if [ $? -ne 0 ]
then
    echo failed INSIDE TICLL-chainclean
    exit
fi


echo "checking chain results...."
diff /tmp/chained.clean $refdir/chained.clean >& /dev/null
if [ $? -ne 0 ]
then
    echo "differences in TICLL-chainclean results"
    echo "using diff /tmp/chained.clean $refdir/chained.clean"
    exit
fi

echo OK
