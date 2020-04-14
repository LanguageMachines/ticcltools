TICCLTOOLS
==================

TICCLtools is a collection of programs to process text data files towards fully-automatic lexical corpus post-correction. Together they constitute the bulk of TICCL: Text Induced Corpus-Cleanup. This software is usually invoked by the pipeline system PICCL: https://github.com/LanguageMachines/PICCL ,
consult there for installation and usage instructions unless you really want to invoke the individual tools manually.

The workflows in PICCL, the Philosophical Integrator of Computational and Corpus Libraries are schematically visualised here, TICCL being the one to the right: 

![PICCL Architecture](https://raw.githubusercontent.com/LanguageMachines/PICCL/master/architecture.png)

Preparation for a specific language and its alphabet:

(Note: A fairly wide range of language specific alphabet and character confusion files are available online, precluding the need for performing this preparatory step yourself.)

- TICCL-lexstat
  - Creates a character frequency ranked 'alphabet' file of the unicode characters that are present in a lexicon for 
    the specific language.
  - Convert an 'alphabet' file (in a second step) into a list of character confusion hash values and an example of the
     particular character confusion, optionally a list of all possible character confusions given the set of characters          involved.
   
   Note that each extra character allowed to be an actual character used in the language expands the search space for           lexical variants. The tool therefore allows you to 'clip' or apply a frequency cut-off to the character frequency 
    list for your particular language.

The actual TICCL post-correction programs in this collection are:
- TICCL-stats
  - A tool to derive word frequency lists from text files or corpora. Its companion tool for corpora in FoLiA XML format,
    FoLiA-stats, is far more developed and recommended.
- TICCL-unk
  - a cleanup tool for word frequency lists. Creates a 'clean' file with desirable words, an 'unk' file with 
    uncorrectable words and a 'punct' file with words that would be clean after removing punctuation before and after.
- TICCL-anahash
  - a tool to create anagram hash values from a word frequency file. All anagrams formable given a particular bag 
     of characters and observed in the frequency file are assigned a distinguishing anagram value, based on the 
     individual character values assigned by tool TICCL-lexstat to each character in the alphabet.
- TICCL-indexer and TICCL-indexerNT
  - a  tool to create an exhaustive numerical index to all lexical
    variants present in a corpus within the distances defined by the character 
    confusion values, given a particular Levenshtein or edit distance.
- TICCL-LDcalc
  - a preprocessing tool for TICCL-rank. Gathers the info from TICCL-anahash, TICCL-indexer or TICCL-indexerNT, 
    TICCL-lexstat and TICCL-unk. Retrieves and pre-filters the symbolic pairs of word variants linked to their 
    Correction Candidates.
- TICCL-rank
  - ranks a word variant list on the basis of a wide range of criteria and the actual set of ranking features 
    specified to be used in the Correction Candidate ranking.
- TICCL-chain
  - Tool designed to gather variants that lie outside the edit distance to the perceived best Correction Candidate, 
    given the Levenshtein distance set earlier in the work flow. This distance is usually two characters.
  - Chaining = "my friends' friends are my friends" : the highest frequency Correction Candidate in the TICCL-rank output 
    list with best-first ranked variants within the set Levenshtein Distance (LD) that act as CCs for further variants 
    beyond this LD (and so on for even greater LDs) is directly linked to these larger LD variants.
- TICCL-chainclean
  - After correction of n-gram frequency lists (with n > 1) word strings may have been assigned different Correction 
    Candidates on the unigram, bi- or trigram correction levels, leading to inconsistencies. This experimental tool tries
    to solve the inconsistencies.
  
  Post-TICCLtools: actual text editing:
  
    We currently only provide for post-editing of texts based on the list of correction candidates collected 
    by TICCLtools for texts or corpora in FoLiA XML. Please see the FoLiA-tools collection for the tool: 
    FoLiA-correct.
