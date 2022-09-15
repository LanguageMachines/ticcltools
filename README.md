[![Project Status: Unsupported – The project has reached a stable, usable state but the author(s) have ceased all work on it. A new maintainer may be desired.](https://www.repostatus.org/badges/latest/unsupported.svg)](https://www.repostatus.org/#unsupported)

# TICCLTOOLS

TICCLtools is a collection of programs to process text data files towards fully-automatic lexical corpus post-correction. Together they constitute the bulk of TICCL: Text Induced Corpus-Cleanup. This software is usually invoked by the pipeline system PICCL: https://github.com/LanguageMachines/PICCL ,
consult there for installation and usage instructions unless you really want to invoke the individual tools manually.

The workflows in PICCL, the Philosophical Integrator of Computational and Corpus Libraries are schematically visualised here, TICCL being the one to the right:

![PICCL Architecture](https://raw.githubusercontent.com/LanguageMachines/PICCL/master/architecture.png)

Preparation for a specific language and its alphabet:

Note: A fairly wide range of language specific alphabet and character confusion files are available online, precluding the need for performing this preparatory step yourself.

We have prepared TICCL for work in many languages, mainly on the basis of available open source lexicons due to Aspell. The language specific files are available here:

* [All languages](http://ticclops.uvt.nl/TICCL.languagefiles.ALLavailable.20160421.tar.gz)
* [Dutch](http://ticclops.uvt.nl/TICCL.languagefiles.nld.20160421.tar.gz)
* [English](http://ticclops.uvt.nl/TICCL.languagefiles.eng.20160421.tar.gz)
* [Finnish](http://ticclops.uvt.nl/TICCL.languagefiles.fin.20160421.tar.gz)
* [French](http://ticclops.uvt.nl/TICCL.languagefiles.fra.20160421.tar.gz)
* [Frisian](http://ticclops.uvt.nl/TICCL.languagefiles.fry.20160421.tar.gz)
* [German](http://ticclops.uvt.nl/TICCL.languagefiles.deu.20160421.tar.gz)
* [German Fraktur](http://ticclops.uvt.nl/TICCL.languagefiles.deu-frak.20160421.tar.gz)
* [Greek (antique)](http://ticclops.uvt.nl/TICCL.languagefiles.grc.20160421.tar.gz)
* [Greek (modern)](http://ticclops.uvt.nl/TICCL.languagefiles.ell.20160421.tar.gz)
* [Icelandic](http://ticclops.uvt.nl/TICCL.languagefiles.isl.20160421.tar.gz)
* [Italian](http://ticclops.uvt.nl/TICCL.languagefiles.ita.20160421.tar.gz)
* [Latin](http://ticclops.uvt.nl/TICCL.languagefiles.lat.20160421.tar.gz)
* [Polish](http://ticclops.uvt.nl/TICCL.languagefiles.pol.20160421.tar.gz)
* [Portuguese](http://ticclops.uvt.nl/TICCL.languagefiles.por.20160421.tar.gz)
* [Romanian](http://ticclops.uvt.nl/TICCL.languagefiles.ron.20160421.tar.gz)
* [Russian](http://ticclops.uvt.nl/TICCL.languagefiles.rus.20160421.tar.gz)
* [Spanish](http://ticclops.uvt.nl/TICCL.languagefiles.spa.20160421.tar.gz)
* [Swedish](http://ticclops.uvt.nl/TICCL.languagefiles.swe.20160421.tar.gz)

Unpack in your main TICCL directory. A subdirectory ``data/int/`` will be
created to house the required files for the specific language(s).

Should you want or need to build your own TICCL alphabet and character confusion files yourself, the tool to do that is:

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
    by TICCLtools for texts or corpora in FoLiA XML. Please see the FoLiA-utils
    (https://github.com/LanguageMachines/foliautils) collection for the tool:
    FoLiA-correct.

## Manual Installation

We provide containers for simple installation, see the next section. If you want to build and install manually on a Linux/BSD system instead, follow these instructions:

First ensure the following dependencies are installed on your system:

* Git
* A sane build environment with a C++ compiler (gcc or clang), make, libtool, pkg-config, autoconf, automake and autoconf-archive
* libbz2, libtar, libicu, libxml2 (including -dev versions for the headers, note that the naming of the packages may vary based on your distribution)
    * On debian/ubuntu, the following should suffice to install the necessary global dependencies: ``sudo apt install make gcc g++ autoconf automake autoconf-archive libtool autotools-dev libicu-dev libxml2-dev libbz2-dev zlib1g-dev libtar-dev``

First ``git clone`` this repository, enter its directory and build as follows:

```console
$ sudo ./build-deps.sh && ./bootstrap.sh && ./configure && make && sudo make install
```

If you have no root permissions, set environment variable ``PREFIX`` to the
target directory where you want to install (ensure it exists), the one in the following example is
a sane default:

```console
$ export PREFIX="$HOME/.local/"
$ ./build-deps.sh && ./bootstrap.sh && ./configure --prefix "$PREFIX" && make && make install
```

Adjust your environment accordingly so the binary and libraries in ``$PREFIX``
can be found: On Linux, ensure the value of ``$PREFIX/lib`` is added to your
`$LD_LIBRARY_PATH` and ``$PREFIX/bin`` directory to your ``$PATH``.

## Container Usage

A pre-made container image can be obtained from Docker Hub as follows:

``docker pull proycon/ticcltools``

You can build a docker container as follows, make sure you are in the root of this repository:

``docker build -t proycon/ticcltools .``

This builds the latest stable release, if you want to use the latest development version
from the git repository instead, do:

``docker build -t proycon/ticcltools --build-arg VERSION=development .``

Run the container interactively as follows:

``docker run -t -i proycon/ticcltools``

Or invoke the tool you want:

``docker run proycon/ticcltools TICCL-rank``

Add the ``-v /path/to/your/data:/data`` parameter (before `-t`) if you want to mount your data volume into the container at `/data`.
