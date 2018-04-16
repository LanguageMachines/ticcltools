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

outdir=OUT/$outsub/TICCL
refdir=OUTreference/BOOK
datadir=DATA
foliadir=BOOK

echo start TICLL-stuff

$bindir/TICCL-lexstat --separator=_ --clip=20 --LD=2 $datadir/nld.aspell.dict

echo "checking lexstat results...."
diff $datadir/nld.aspell.dict.clip20.lc.chars $refdir/dict.lc.chars > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-lexstat .lc results"
    echo "using: diff $datadir/nld.aspell.dict.clip20.lc.chars $refdir/dict.lc.chars"
    exit
fi

diff $datadir/nld.aspell.dict.clip20.ld2.charconfus $refdir/charconfus > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-lexstat .confusion results"
    echo "using: diff $datadir/nld.aspell.dict.clip20.ld2.charconfus $refdir/charconfus"
    exit
fi

echo "start FoLiA-stats..."

$bindir/FoLiA-stats -R -s -t 30 -e folia.xml$ --lang=none --class=OCR --ngram 1 -o $outdir/TESTDP035 --hemp=$outdir/TESTDP035.hemp $foliadir

if [ $? -ne 0 ]
then
    echo failed after FoLiA-stats
    exit
fi

echo start TICLL-unk

cp $outdir/TESTDP035.wordfreqlist.tsv $outdir/TESTDP035.tsv

$bindir/TICCL-unk --background $datadir/INLandAspell.corpus --artifrq 100000000 --acro $outdir/TESTDP035.tsv

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
    echo "using: diff $outdir/TESTDP035.tsv.punct $refdir/punct"
    exit
fi
diff $outdir/TESTDP035.tsv.unk $refdir/unk > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK unk results"
    echo "using: diff $outdir/TESTDP035.tsv.unk $refdir/unk"
    exit
fi
diff $outdir/TESTDP035.tsv.clean $refdir/clean > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK clean results"
    echo "using: diff $outdir/TESTDP035.tsv.clean $refdir/clean"
    exit
fi
diff $outdir/TESTDP035.tsv.acro $refdir/acro > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "differences in Ticcl-UNK acro results"
    echo "using: diff $outdir/TESTDP035.tsv.acro $refdir/acro"
    exit
fi

echo "start TICLL-anahash"

$bindir/TICCL-anahash --alph $datadir/nld.aspell.dict.clip20.lc.chars --artifrq 100000000 $outdir/TESTDP035.tsv.clean

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
    echo "using: diff /tmp/foci $refdir/foci"
    exit
fi

echo "start TICCL-indexer"

$bindir/TICCL-indexer -t 6 --hash $outdir/TESTDP035.tsv.clean.anahash --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus --foci $outdir/TESTDP035.tsv.clean.corpusfoci

if [ $? -ne 0 ]
then
    echo "failed in TICCL-indexer"
    exit
fi

echo "checking INDEXER results...."
diff $outdir/TESTDP035.tsv.clean.index $refdir/index > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-indexer results"
    echo "using diff $outdir/TESTDP035.tsv.clean.index $refdir/index"
    exit
fi

echo "start TICLL-LDcalc"

$bindir/TICCL-LDcalc --index $outdir/TESTDP035.tsv.clean.index --hash $outdir/TESTDP035.tsv.clean.anahash --clean $outdir/TESTDP035.tsv.clean --LD 2 -t 1 --artifrq 100000000 -o $outdir/TESTDP035.tsv.clean.ldcalc

if [ $? -ne 0 ]
then
    echo "failed in TICCL-LDcalc"
    exit
fi

echo "checking LDCALC results...."
sort $outdir/TESTDP035.tsv.clean.ldcalc  > /tmp/ldcalc

diff /tmp/ldcalc $refdir/ldcalc > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "differences in Ticcl-ldcalc results"
    echo "using: diff /tmp/ldcalc $refdir/ldcalc"
    exit
fi

echo "start TICLL-rank"

$bindir/TICCL-rank -t 30 --alph $datadir/nld.aspell.dict.clip20.lc.chars --charconf $datadir/nld.aspell.dict.clip20.ld2.charconfus -o $outdir/TESTDP035.tsv.clean.ldcalc.ranked --debugfile $outdir/.TESTDP035.tsv.clean.ldcalc.debug.ranked --artifrq 0 --clip 5 --skipcols=10,11 $outdir/TESTDP035.tsv.clean.ldcalc 2> $outdir/.TESTDP035.RANK.stderr

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
    echo "using: diff /tmp/rank.sorted $refdir/rank.sorted"
    exit
else
    echo OK!
fi
