Run TICCLopsTEST.pl with the command line:

perl /reddata/TICCLAT/TEST/TICCLopsTEST.pl /reddata/TICCLAT/TEST/TICCLconfig/TICCL.TEST.config

This will run and test the full TICCL suite of C++ modules.

If the output directory OUT/ exists, files in it will be overwritten. You may either remove it first, or specify another ouput directory in configuration file /reddata/TICCLAT/TEST/TICCLconfig/TICCL.TEST.config.

The directory OUTreference/ is a copy of the OUT/ directory produced for reference purposes in an earlier test run.

The script now outputs to stderr the actual command line used for each C++ module, for diagnostic purposes. E.g.:

RUN_FoLiA-stats1FOLIA: Command line: $ /exp/sloot/usr/local/bin//FoLiA-stats -R -s -t 30 -e folia.xml$ --lang=none --class=OCR --ngram 1 -o /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035 /reddata/TICCLAT/TEST/FOLIA/
RUN_TICCL-unk: Command line: $ /exp/sloot/usr/local/bin//TICCL-unk --corpus /exp2/reynaert/ticclops/data/int/nld/nuTICCL.OldandINLlexandINLNamesAspell.v2.COL1.tsv --artifrq 100000000 /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv
RUN_TICCL-anahash: Command line: $ /exp/sloot/usr/local/bin//TICCL-anahash --alph /exp2/reynaert/ticclops/data/int/nld/nld.aspell.dict.lc.chars --artifrq 100000000 /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean
RUN_TICCL-indexerNT: Command line: $ /exp/sloot/usr/local/bin//TICCL-indexerNT -t 30 --hash /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.anahash --charconf /exp2/reynaert/ticclops/data/int/nld/nld.aspell.dict.c20.d2.confusion --foci /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.corpusfoci -o .tsv.clean.confuslist
RUN_TICCL-LDcalc-withm: Command line: $ /exp/sloot/usr/local/bin//TICCL-LDcalc --index .tsv.clean.confuslist.indexNT --hash /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.anahash --clean /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean --LD 2 -t 30 --artifrq 100000000 -o /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.ldcalc
RUN_TICCL-rank2: Command line: $ /exp/sloot/usr/local/bin//TICCL-rank -t 30 --alph /exp2/reynaert/ticclops/data/int/nld/nld.aspell.dict.lc.chars --charconf /exp2/reynaert/ticclops/data/int/nld/nld.aspell.dict.c20.d2.confusion -o /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.ldcalc.ranked --debugfile /reddata/TICCLAT/TEST/OUT//zzz/TICCL/.TESTDP035.tsv.clean.ldcalc.debug.ranked --artifrq 0 --clip 5 --skipcols=10,11 /reddata/TICCLAT/TEST/OUT//zzz/TICCL/TESTDP035.tsv.clean.ldcalc 2>/reddata/TICCLAT/TEST/OUT//zzz/TICCL/.TESTDP035.RANK.stderr

MRE
2017-11-20
