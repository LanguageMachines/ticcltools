# $Id$
# $URL$

AM_CPPFLAGS = -I@top_srcdir@/include
AM_CXXFLAGS = -std=c++11 -g -W -Wall -pedantic -O3

bin_PROGRAMS = TICCL-indexer TICCL-indexerNT \
	TICCL-LDcalc TICCL-unk TICCL-lexstat \
	TICCL-anahash TICCL-rank TICCL-lexclean \
	W2V-near W2V-dist

LDADD = libticcl.la
lib_LTLIBRARIES = libticcl.la
libticcl_la_LDFLAGS= -version-info 1:0:0

libticcl_la_SOURCES = word2vec.cxx

TICCL_indexer_SOURCES = TICCL-indexer.cxx
TICCL_indexerNT_SOURCES = TICCL-indexerNT.cxx
#indexerNT_SOURCES = indexerNT.cxx
#indexerNT2_SOURCES = indexerNT2.cxx
TICCL_LDcalc_SOURCES = TICCL-LDcalc.cxx
TICCL_rank_SOURCES = TICCL-rank.cxx
TICCL_unk_SOURCES = TICCL-unk.cxx
TICCL_lexstat_SOURCES = TICCL-lexstat.cxx
TICCL_lexclean_SOURCES = TICCL-lexclean.cxx
TICCL_anahash_SOURCES = TICCL-anahash.cxx
W2V_near_SOURCES = W2V-near.cxx
W2V_dist_SOURCES = W2V-dist.cxx
