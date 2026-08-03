#!/usr/bin/env ruby
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
# match with an E value of less than 1e−5, a sequence identity of more than 40%, and an alignable region
# of the query protein sequence of more than 50%. For a pair of genomes, each genome was used as the
# query genome to perform the BLASTP search. The number of conserved proteins in each genome of strains
# being compared was slightly different because of the existence of duplicate genes (paralogs). The percentage
# of conserved proteins (POCP) between two genomes was calculated as [(C1 + C2)/(T1 + T2)] · 100%, where C1 and
# C2 represent the conserved number of proteins in the two genomes being compared, respectively, and T1 and T2
# represent the total number of proteins in the two genomes being compared, respectively. In theory, the POCP
# value can vary from 0% to 100%, depending on the similarity of the protein contents of two genomes.

########################################################################

## produce a matrix/excel format from the POCP pairwise output files
strains = {}

Dir.glob("*-vs-*.txt").each do |comp|
  bn = File.basename(comp, '.txt')
  g1, g2 = bn.split('-vs-')
  strains[g1] = [] unless strains.key?(g1)
  strains[g2] = [] unless strains.key?(g2)
end
puts "Collected #{strains.keys.size} strains."

## sort the IDs so that the matrix layout is the same for every run
ids = strains.keys.sort

ids.each do |strain1|
  ids.each do |strain2|
    if strain1 == strain2
      strains[strain1].push('100.0')
    else
      comp = ["#{strain1}-vs-#{strain2}.txt", "#{strain2}-vs-#{strain1}.txt"].find { |f| File.exist?(f) }
      if comp
        strains[strain1].push(File.read(comp).chomp.strip.to_f.round(4).to_s)
      else
        ## in one-vs-all mode most pairs are not calculated at all, report them as
        ## missing instead of silently writing a 0.0 POCP value
        strains[strain1].push('NA')
      end
    end
  end
end

out = File.open("pocp-matrix.tsv",'w')
out << "ID\t" << ids.join("\t") << "\n"
ids.each do |strain_id|
  out << "#{strain_id}\t" << strains[strain_id].join("\t") << "\n"
end
out.close