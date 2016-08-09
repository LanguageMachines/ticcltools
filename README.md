TICCLTOOLS

TicclTools is a colection of programs to process datafiles as used in Ticcl.

These datafiles are in general large, ans sometimes complex.

The main programs in thsi colection are:
- TICCL-indexer and TICCL-indexerNT
  - a  tool  to create an exhaustive index to all lexical
    variants given a particular Levenshtein or edit distance in a corpus.
- TICCL-anahash
  - a tool to create anagram hashes form a word frequency file. Also creates ab 'alphabet' file of the unicode characters that are present in the corpus.
- TICLL-LDcalc
  - a proprocessing tool for TICCL-rank. Gathers the info from TICC-anahash, TICCL-indexer, TICCL-lexstat and TICCL-unk
- TICCL-rank
  - ranks a word varian list according to al lot of criteria
- TICCL-unk
  - a cleanup tool for word frequency lists. creates a 'clean' file with desirable words, an 'unk' file with uncorrectable words and a 'punct' file with words that would be clean after removing puncuation.
- TICLL-lexstat
  - convert an 'alphabet' file (from TICCL-anahash) into a frequency list of hashes and optionally a list of confusions.

