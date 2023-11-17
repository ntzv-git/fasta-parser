# fasta-parser

<pre>
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8+
Version : 1.0
</pre>


## Description

Analyzes fasta36 tools with output format '-m 0' and returns a tabular similar to output format '-m 8', adding the
aligned sequences of query (q_aln) and subject (s_aln) as well as the matched aligned pattern (m_aln).
The matched aligned pattern is a string composed by either ":" for a perfect match, "." when the mismatch score >= 0
or " " (white space) for other cases.

Inspired by "parse_ssearch.py" from https://github.com/jtremblay/MiRNATarget


## Installation

```bash
git clone https://github.com/ntzv-git/fasta-parser.git
```


### Requirements

- python3.8


### Package dependencies

- argparse
- re
- sys


## Usage

```
usage: fasta-parser.py [-h] -i INPUT -o OUTPUT

fasta-parser.py v1.0 parses the "-m 0" output format of the fasta36 tool and
returns an table similar to the "-m 8" output format, including the aligned
pattern and sequences

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input alignment file
  -o OUTPUT, --output OUTPUT
                        path to output tabular file
```
