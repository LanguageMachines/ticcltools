#!/bin/bash

if [ "$1" != "" ]
then
    outsub=$1
else
    outsub=zzz
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

$bindir/TICCL-lexstat --clip=20 --LD=2 -o $outdir/aspell $datadir/nld.aspell.dict

echo "checking lexstat results...."
diff $outdir/aspell.lc.chars $refdir/dict.lc.chars > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-lexstat .lc results"
    echo "using: diff $outdir/aspell.lc.chars $refdir/dict.lc.chars"
    exit
fi

diff $outdir/aspell.clip20.ld2.charconfus $refdir/charconfus > /dev/null 2>&1
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

$bindir/TICCL-unk --corpus $datadir/INLandAspell.corpus --artifrq 100000000 -o $outdir/TESTDP035 $outdir/TESTDP035.wordfreqlist.tsv

if [ $? -ne 0 ]
then
    echo "failed in TICCL-unk"
    exit
fi
echo "checking UNK results...."
diff $outdir/TESTDP035.punct $refdir/punct > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK punct results"
    echo "using: diff $outdir/TESTDP035.wordfreqlist.tsv.punct $refdir/punct"
    exit
fi

echo "start TICLL-anahash"

$bindir/TICCL-anahash --alph $outdir/aspell.lc.chars --artifrq 100000000 $outdir/TESTDP035.clean

if [ $? -ne 0 ]
then
    echo failed after TICCL-anahash
    exit
fi

echo "checking ANAHASH results...."
sort $outdir/TESTDP035.clean.corpusfoci > /tmp/foci
diff /tmp/foci $refdir/foci > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-anahash foci results"
    exit
fi

echo "start TICCL-indexerNT"

$bindir/TICCL-indexerNT -t 30 --hash $outdir/TESTDP035.clean.anahash --charconf $outdir/aspell.clip20.ld2.charconfus --foci $outdir/TESTDP035.clean.corpusfoci

if [ $? -ne 0 ]
then
    echo "failed in TICCL-indexerNT"
    exit
fi

echo "checking INDEXER results...."
sort $outdir/TESTDP035.clean.indexNT > /tmp/indexNT
diff /tmp/indexNT $refdir/indexNT > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-indexer results"
    echo "using diff /tmp/indexNT $refdir/indexNT"
    exit
fi

echo "start TICLL-LDcalc"

$bindir/TICCL-LDcalc --index $outdir/TESTDP035.clean.indexNT --hash $outdir/TESTDP035.clean.anahash --clean $outdir/TESTDP035.clean --LD 2 -t 30 --artifrq 100000000 -o $outdir/TESTDP035.clean.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICCL-LDcalc"
    exit
fi

echo "checking LDCALC results...."
sort $outdir/TESTDP035.clean.ldcalc  > /tmp/ldcalc

diff /tmp/ldcalc $refdir/ldcalc > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    ech "using: diff /tmp/ldcalc $refdir/ldcalc"
    exit
fi

echo "start TICLL-rank"

$bindir/TICCL-rank -t 30 --alph $outdir/aspell.lc.chars --charconf $outdir/aspell.clip20.ld2.charconfus -o $outdir/TESTDP035.ldcalc.ranked --debugfile $outdir/.TESTDP035.ldcalc.debug.ranked --artifrq 0 --clip 5 --skipcols=10,11 $outdir/TESTDP035.clean.ldcalc 2> $outdir/.TESTDP035.RANK.stderr

if [ $? -ne 0 ]
then
    echo failed after TICLL-rank
    exit
fi

echo "checking RANK results...."

sort $outdir/TESTDP035.ldcalc.ranked > /tmp/rank.sorted
diff /tmp/rank.sorted $refdir/rank.sorted > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in TICLL-rank results"
    exit
else
    echo OK!
fi
