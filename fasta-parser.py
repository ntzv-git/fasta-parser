#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8
Date    : 18/11/2023

Description :
Analyzes fasta36 tools with output format '-m 0' and returns a tabular similar to output format '-m 8', adding the
aligned sequences of query (q_aln) and subject (s_aln) as well as the matched aligned pattern (m_aln).
The matched aligned pattern is a string composed by either ":" for a perfect match, "." when the mismatch score >= 0
or " " (white space) for other cases.
Inspired by "parse_ssearch.py" from https://github.com/jtremblay/MiRNATarget

LICENSE :
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import argparse
import re
import sys

VERSION = "1.0"


def parse_arguments():
    """
    Parse the system arguments.
    """
    parser = argparse.ArgumentParser(description=f'fasta-parser.py v{VERSION} parses the "-m 0" output format of the '
                                                 f'fasta36 tool and returns an table similar to the "-m 8" output '
                                                 f'format, including the aligned pattern and sequences')
    parser.add_argument('-i', '--input', required=True, help='path to input alignment file')
    parser.add_argument('-o', '--output', required=True, help='path to output tabular file')
    args = parser.parse_args()
    return args


def count_mutations(seq1: str, seq2: str) -> (int, int, int):
    """
    Count mutations between seq1 and seq2 (they must have the same length).
    :param seq1: sequence 1
    :param seq2: sequence 3
    :return: Numbre of mismatches, gaps and opening gaps between seq1 and seq2
    """
    gaps = seq1.count("-") + seq2.count("-")

    gap_opens = len(re.findall(r"[\w\*]-", seq1)) + len(re.findall(r"[\w\*]-", seq2))
    if seq1[0] == "-":
        gap_opens += 1
    if seq2[0] == "-":
        gap_opens += 1

    mismatches = 0
    for i, char in enumerate(seq1):
        # Do not use gaps to count mismatches
        if (seq1[i] != "-") and (seq2[i] != "-") and (seq1[i] != seq2[i]):
            mismatches += 1

    return mismatches, gaps, gap_opens


def write(outfile, cursor: int, query: str, subject: str, p_ident: str, aln_len: str, q_start: str, q_end: str,
          s_start: str, s_end: str, e_value: str, bit_score: str, p_sim: str, q_len: str, s_len: str,
          m_aln: str, q_aln: str, s_aln: str) -> None:
    """
    Write the output line only if all data have been recovered.
    """
    if cursor == 3:
        mismatches, gaps, gap_opens = count_mutations(q_aln, s_aln)
        outfile.write(f"{query}\t{subject}\t{p_ident}\t{aln_len}\t{mismatches}\t{gap_opens}\t{q_start}\t{q_end}\t"
                      f"{s_start}\t{s_end}\t{e_value}\t{bit_score}\t{p_sim}\t{gaps}\t{q_len}\t{s_len}\t"
                      f"{m_aln}\t{q_aln}\t{s_aln}\n")


def main():
    # Import parameters
    args = parse_arguments()
    inpath = args.input
    outpath = args.output

    print(f"inpath  : {inpath}")
    print(f"outpath : {outpath}")
    print("")

    # Load input file
    infile = open(inpath, 'rt')
    inlines = infile.readlines()
    infile.close()

    # Initilaization
    cursor = 0
    query = subject = q_len = s_len = bit_score = e_value = p_ident = p_sim = aln_len = q_start = q_end = \
        s_start = s_end = q_aln = s_aln = m_aln = ""
    match_span = (-1, -1)
    parse_pattern = False
    outfile = open(outpath, 'wt')
    outfile.write("#query\tsubject\tp_ident\taln_len\tmismatches\tgap_opens\tq_start\tq_end\ts_start\ts_end\tevalue\t"
                  "bit_score\tp_sim\tgaps\tq_len\ts_len\tm_aln\tq_aln\ts_aln\n")

    # Processing
    for line in inlines:
        line = line[:-1]

        # print(f"line : {line}")
        # if line == ">>XM_048752055.1_frame1                                   (1176 aa)":
        #     break

        # Query title
        match = re.match(r"^\s*\d+>>>(\S+) - (\d+) [na][ta]$", line)
        if match:
            write(outfile, cursor, query, subject, p_ident, aln_len, q_start, q_end, s_start, s_end, e_value,
                  bit_score, p_sim, q_len, s_len, m_aln, q_aln, s_aln)
            cursor = 1
            query = match.group(1)
            q_len = match.group(2)
            subject = s_len = ""
            continue

        # Subject title
        match = re.match(r"^>>(\S+) .*\((\d+) [na][ta]\)$", line)
        if match:
            write(outfile, cursor, query, subject, p_ident, aln_len, q_start, q_end, s_start, s_end, e_value,
                  bit_score, p_sim, q_len, s_len, m_aln, q_aln, s_aln)
            cursor = 2
            subject = match.group(1)
            s_len = match.group(2)
            bit_score = e_value = p_ident = p_sim = aln_len = q_start = q_end = s_start = s_end = \
                q_aln = s_aln = m_aln = ""
            match_span = (-1, -1)
            continue

        # When the process concerns the same subject (equivalent to the subject title)
        if (line == ">--") and (cursor == 3):
            write(outfile, cursor, query, subject, p_ident, aln_len, q_start, q_end, s_start, s_end, e_value,
                  bit_score, p_sim, q_len, s_len, m_aln, q_aln, s_aln)
            cursor = 2
            bit_score = e_value = p_ident = p_sim = aln_len = q_start = q_end = s_start = s_end = \
                q_aln = s_aln = m_aln = ""
            match_span = (-1, -1)
            continue

        # Alignment score
        match = re.match(r"^ .+ bits: (.+) E\(\d+\): (.+)$", line)
        if match and (cursor == 2):
            bit_score = match.group(1)
            e_value = match.group(2)
            continue

        # Alignment details
        match = re.match(r"^\w.+score: \d+; (.+)% identity \((.+)% similar\) in (\d+) [na][ta] overlap "
                         r"\((\d+)-(\d+):(\d+)-(\d+)\)$", line)
        if match and (cursor == 2):
            cursor = 3
            p_ident = match.group(1)
            p_sim = match.group(2)
            aln_len = match.group(3)
            q_start = match.group(4)
            q_end = match.group(5)
            s_start = match.group(6)
            s_end = match.group(7)
            continue

        # Alignment graph
        if cursor == 3:
            # if the line start by the query name
            if line[:4] == query[:4]:
                match = re.match(r"^\S+\s+(\S+)\s*$", line)
                match_span = match.span(1)
                q_aln += line[match_span[0]:match_span[1]]
                parse_pattern = True
            # if the line start by the subject name
            elif line[:4] == subject[:4] and parse_pattern:
                s_aln += line[match_span[0]:match_span[1]]
                parse_pattern = False
            # parse next line after query
            elif parse_pattern:
                m_aln += line[match_span[0]:match_span[1]]

    write(outfile, cursor, query, subject, p_ident, aln_len, q_start, q_end, s_start, s_end, e_value,
          bit_score, p_sim, q_len, s_len, m_aln, q_aln, s_aln)
    outfile.close()


if __name__ == '__main__':
    main()
