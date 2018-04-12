#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=roar
fi

bindir=/home/sloot/usr/local/bin

if [ !-e $bindir ]
then
   bindir=/exp/sloot/usr/local/bin
   if [ !-e $bindir ]
   then
       echo "cannot find executables "
       exit
   fi
fi
outdir=OUT/$outsub/TICCL
refdir=OUTreference/OK
datadir=DATA

echo start TICLL-stuff

$bindir/TICCL-lexstat --clip=20 --LD=2 $datadir/nld.aspell.dict

echo "checking lexstat results...."
diff $datadir/nld.aspell.dict.lc.chars $refdir/dict.lc.chars > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-lexstat .lc results"
    exit
fi

diff $datadir/nld.aspell.dict.clip20.ld2.charconfus $refdir/charconfus > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-lexstat .confusion results"
    exit
fi

echo "start FoLiA-stats..."

$bindir/FoLiA-stats -R -s -t 30 -e folia.xml$ --lang=none --class=OCR --ngram 1 -o $outdir/TESTDP035 --hemp=$outdir/TESTDP035.hemp FOLIA/

if [ $? -ne 0 ]
then
    echo failed after FoLiA-stats
    exit
fi

echo start TICLL-unk

cp $outdir/TESTDP035.wordfreqlist.tsv $outdir/TESTDP035.tsv

$bindir/TICCL-unk --background $datadir/INLandAspell.corpus --artifrq 100000000 $outdir/TESTDP035.tsv

if [ $? -ne 0 ]
then
    echo "failed in TICCL-unk"
    exit
fi
echo "checking UNK results...."
diff $outdir/TESTDP035.tsv.punct $refdir/punct > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK punct results"
    exit
fi

echo "start TICLL-anahash"

$bindir/TICCL-anahash --alph $datadir/nld.aspell.dict.lc.chars --artifrq 100000000 $outdir/TESTDP035.tsv.clean

if [ $? -ne 0 ]
then
    echo failed after TICCL-anahash
    exit
fi

echo "checking ANAHASH results...."
sort $outdir/TESTDP035.tsv.clean.corpusfoci > /tmp/foci
diff /tmp/foci $refdir/foci > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-anahash foci results"
    exit
fi

echo "start TICCL-indexerNT-roaring"

$bindir/TICCL-indexerNT-roaring -t 30 --hash $outdir/TESTDP035.tsv.clean.anahash --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus --foci $outdir/TESTDP035.tsv.clean.corpusfoci

if [ $? -ne 0 ]
then
    echo "failed in TICCL-indexerNT-roaring"
    exit
fi

echo "checking INDEXER results...."
sort $outdir/TESTDP035.tsv.clean.indexNT.R > /tmp/indexNT.R
diff /tmp/indexNT.R $refdir/indexNT.R > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-indexer results"
    echo "using diff /tmp/indexNT $refdir/indexNT"
    exit
fi

echo "start TICLL-LDcalc-roaring"

$bindir/TICCL-LDcalc-roaring --index $outdir/TESTDP035.tsv.clean.indexNT.R --hash $outdir/TESTDP035.tsv.clean.anahash --clean $outdir/TESTDP035.tsv.clean --LD 2 -t 30 --artifrq 100000000 -o $outdir/TESTDP035.tsv.clean.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICCL-LDcalc-roaring"
    exit
fi

echo "checking LDCALC results...."
sort $outdir/TESTDP035.tsv.clean.ldcalc  > /tmp/ldcalc

diff /tmp/ldcalc $refdir/ldcalc > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    ech "using: diff /tmp/ldcalc $refdir/ldcalc"
    exit
fi

echo "start TICLL-rank"

$bindir/TICCL-rank -t 30 --alph $datadir/nld.aspell.dict.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/TESTDP035.tsv.clean.ldcalc.ranked --debugfile $outdir/.TESTDP035.tsv.clean.ldcalc.debug.ranked --artifrq 0 --clip 5 --skipcols=10,11 $outdir/TESTDP035.tsv.clean.ldcalc 2> $outdir/.TESTDP035.RANK.stderr

if [ $? -ne 0 ]
then
    echo failed after TICLL-rank
    exit
fi

echo "checking RANK results...."

sort $outdir/TESTDP035.tsv.clean.ldcalc.ranked > /tmp/rank.sorted
diff /tmp/rank.sorted $refdir/rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    exit
else
    echo OK!
fi
