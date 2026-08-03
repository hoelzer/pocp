#!/usr/bin/env python3
#
# Author: hoelzer.martin@gmail.com

###################################
## P.O.C.P - calculate percentage of conserved proteins
##
## A prokaryotic genus can be defined as a group of species
## with all pairwise POCP values higher than 50%
##
## http://jb.asm.org/content/196/12/2210.full

## POCP. The conserved proteins between a pair of genomes were determined by aligning all the protein
# sequences of one genome (query genome) with all the protein sequences of another genome using the
# BLASTP program. Proteins from the query genome were considered conserved when they had a BLAST
# match with an E value of less than 1e-5, a sequence identity of more than 40%, and an alignable region
# of the query protein sequence of more than 50%. For a pair of genomes, each genome was used as the
# query genome to perform the BLASTP search. The number of conserved proteins in each genome of strains
# being compared was slightly different because of the existence of duplicate genes (paralogs). The percentage
# of conserved proteins (POCP) between two genomes was calculated as [(C1 + C2)/(T1 + T2)] * 100%, where C1 and
# C2 represent the conserved number of proteins in the two genomes being compared, respectively, and T1 and T2
# represent the total number of proteins in the two genomes being compared, respectively. In theory, the POCP
# value can vary from 0% to 100%, depending on the similarity of the protein contents of two genomes.

########################################################################

## produce a matrix/excel format from the POCP pairwise output files

import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Combine the pairwise POCP values into a matrix.')
parser.add_argument('--output', default='pocp-matrix.tsv', help='TSV file to write the matrix to')
args = parser.parse_args()

strains = set()
for comparison in glob.glob("*-vs-*.txt"):
    g1, g2 = os.path.basename(comparison)[:-len('.txt')].split('-vs-')
    strains.add(g1)
    strains.add(g2)
print(f"Collected {len(strains)} strains.")

## sort the IDs so that the matrix layout is the same for every run
ids = sorted(strains)

matrix = {}
for strain1 in ids:
    row = []
    for strain2 in ids:
        if strain1 == strain2:
            row.append('100.0')
            continue

        comparison = next(
            (f for f in (f"{strain1}-vs-{strain2}.txt", f"{strain2}-vs-{strain1}.txt") if os.path.exists(f)),
            None
        )
        if comparison:
            with open(comparison) as fh:
                row.append(str(round(float(fh.read().strip()), 4)))
        else:
            ## in one-vs-all mode most pairs are not calculated at all, report them as
            ## missing instead of silently writing a 0.0 POCP value
            row.append('NA')
    matrix[strain1] = row

with open(args.output, 'w') as out:
    out.write("ID\t" + "\t".join(ids) + "\n")
    for strain_id in ids:
        out.write(f"{strain_id}\t" + "\t".join(matrix[strain_id]) + "\n")
