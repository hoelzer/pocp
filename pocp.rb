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


class Pocp

  EVALUE = 1e-5
  SEQ_IDENTITY = 0.4
  ALN_LENGTH = 0.5

  def initialize(species1_faa, species2_faa, out_dir, threads)

    puts "##################################\n## P.O.C.P Blaster \n##################################\n\nStart:\t#{Time.now}\n\n"
    puts "Species 1:\t#{File.basename(species2_faa)}"
    puts "Species 2:\t#{File.basename(species1_faa)}"
    puts "Write output: #{out_dir}\n\n"

    species1_bn = File.basename(species1_faa)
    species2_bn = File.basename(species2_faa)

    Dir.mkdir(out_dir) unless Dir.exists?(out_dir)
    blast_dir = "#{out_dir}/blast_db/"
    Dir.mkdir(blast_dir) unless Dir.exists?(blast_dir)
    `cp #{species2_faa} #{blast_dir}`
    `cp #{species1_faa} #{blast_dir}`
    species2_faa = "#{blast_dir}/#{File.basename(species2_faa)}"
    species1_faa = "#{blast_dir}/#{File.basename(species1_faa)}"

    # build blast target db unless exists
    puts 'Make blast databases...'
    [species1_faa, species2_faa].each do |faa|
      if File.exists?("#{faa}.psi")
        puts "\tAlready exists for #{File.basename(faa)}!"
      else
        `makeblastdb -in #{faa} -dbtype prot -parse_seqids`
        puts "\tDone for #{File.basename(faa)}."
      end
    end

    # blast
    # run blast only once per genome pair
    puts "Run Blastp..."
    c1 = blast(species1_faa, species2_faa, species1_bn, out_dir, threads)
    c2 = blast(species2_faa, species1_faa, species2_bn, out_dir, threads)
    t1 = `grep ">" #{species1_faa} | wc -l`.to_i
    t2 = `grep ">" #{species2_faa} | wc -l`.to_i

    pocp = (((c1.to_f + c2)/(t1 + t2)) * 100).round(4)

    puts "\n\nC1:\t#{c1}\nC2:\t#{c2}\n\nT1:\t#{t1}\nT2:\t#{t2}\n\nPOCP = [(C1+C2)/(T1+T2)]*100% = #{pocp}"
    puts "\nEnd:\t#{Time.now}"

    # write final output
    f = File.open("#{out_dir}/pocp.txt",'w')
    f << "##################################\n## P.O.C.P Blaster \n##################################\n\nStart:\t#{Time.now}\n\n"
    f << "Species 1:\t#{File.basename(species2_faa)}\n"
    f << "Species 2:\t#{File.basename(species1_faa)}\n"
    f << "\n\nC1:\t#{c1}\nC2:\t#{c2}\n\nT1:\t#{t1}\nT2:\t#{t2}\n\nPOCP = [(C1+C2)/(T1+T2)]*100% = #{pocp}"
    f << "\n\nEnd:\t#{Time.now}"
    f.close
  end

  def blast(species1_faa, species2_faa, species1_bn, out_dir, threads)
    `blastp -task blastp -num_threads #{threads} -query #{species1_faa} -db #{species2_faa} -evalue #{EVALUE} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" | awk '{if($3>#{SEQ_IDENTITY*100} && $4>($9*#{ALN_LENGTH})){print $0}}' > #{out_dir}/#{species1_bn}.blast`
    c = `awk '{print $1}' #{out_dir}/#{species1_bn}.blast | sort | uniq | wc -l`.split(' ')[0].to_i
    puts "\t#{species1_bn}:\tFound #{c} matches with an E value of less than #{EVALUE}, a sequence identity of more than #{SEQ_IDENTITY*100}%, and an alignable region of the query protein sequence of more than #{ALN_LENGTH*100}%."
    c
  end
end


fasta_dir = ARGV[0]
out_dir = ARGV[1]
threads = ARGV[2]

strains = []

Dir.glob(fasta_dir+"/*.faa").each do |fasta1|
  species1_faa = fasta1
  species1_name = File.basename(fasta1, '.faa')
  Dir.glob(fasta_dir+"/*.faa").each do |fasta2|
    species2_faa = fasta2
    species2_name = File.basename(fasta2, '.faa')
    out = "#{out_dir}/#{species1_name}-vs-#{species2_name}"
    `mkdir -p #{out}`
    Pocp.new(species1_faa, species2_faa, out, threads)
    strains.push(species1_name) unless strains.include?(species1_name)
    strains.push(species2_name) unless strains.include?(species2_name)
  end
end

#ruby $DIR/matrix.rb "/home/hoelzer/projects/fabien_rsha_pocp/calc/flamingo" "14-2711_R47,15-2067_O50,15-2067_O99,Cav_10DC88,Cga_08-1274-3,Cps_6BC,Cab_S26-3,Cca_GPIC,Cfe_Fe-C56,Cpe_E58,Cpn_TW-183,Ctr_D-UW-3-CX,Cmu_Nigg,Csu_SWA-2,Cib_10-1398-6"
## produce a matrix/excel format from the POCP script output
pocps = {}

Dir.glob("#{out_dir}/*-vs-*").each do |comp|
	bn = File.basename(comp)
  g1 = bn.split('-vs-')[0]
  g2 = bn.split('-vs-')[1]

  pocp = `grep POCP #{comp}/pocp.txt | awk 'BEGIN{FS=" "};{print $5}'`.chomp.strip

  if pocps[g1]
	  pocps[g1].push(g2,pocp)
  else
	  pocps[g1] = [g2,pocp]
  end  
  strains.push(g1) unless strains.include?(g1)
  strains.push(g2) unless strains.include?(g2)
end
puts "read in POCP for #{strains.size} strains."

out = File.open("#{out_dir}/#{File.basename(out_dir)}.csv",'w')
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
	    comp = "#{strain1}-vs-#{strain2}"
		  pocp = `grep POCP #{out_dir}/#{comp}/pocp.txt | awk 'BEGIN{FS=" "};{print $5}'`.chomp.strip
			out_line += ",#{pocp}"
    end    
  end
  out << out_line << "\n"  
end
out.close

