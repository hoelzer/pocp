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

#ruby $DIR/matrix.rb "/home/hoelzer/projects/fabien_rsha_pocp/calc/flamingo" "14-2711_R47,15-2067_O50,15-2067_O99,Cav_10DC88,Cga_08-1274-3,Cps_6BC,Cab_S26-3,Cca_GPIC,Cfe_Fe-C56,Cpe_E58,Cpn_TW-183,Ctr_D-UW-3-CX,Cmu_Nigg,Csu_SWA-2,Cib_10-1398-6"

## produce a matrix/excel format from the POCP script output
pocps = {}
strains = []

Dir.glob("*-vs-*").each do |comp|
	bn = File.basename(comp)
  g1 = bn.split('-vs-')[0]
  g2 = bn.split('-vs-')[1].sub('.txt','')

  pocp = `cat #{comp}`.chomp.strip.to_f.round(4)

  if pocps[g1]
	  pocps[g1].push(g2,pocp)
  else
	  pocps[g1] = [g2,pocp]
  end  
  strains.push(g1) unless strains.include?(g1)
  strains.push(g2) unless strains.include?(g2)
end
puts "read in POCP for #{strains.size} strains."

out = File.open("pocp-matrix.csv",'w')
out << 'ID,' << strains.join(',') << "\n"

line_count = 0
strains.each do |strain1|
  out_line = strain1
  line_count += 1
  tmp_count = line_count
	strains.each do |strain2|
    tmp_count = tmp_count - 1
    if tmp_count > 0 
		  out_line += ","
    else      
      comp = "#{strain1}-vs-#{strain2}.txt"
      if File.exist?(comp)
		    pocp = `cat #{comp}`.chomp.strip.to_f.round(4)
        out_line += ",#{pocp}"
      end
    end    
  end
  out << out_line << "\n"  
end
out.close