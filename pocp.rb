#!/home/hoelzer/local/bin/ruby

###################################
## P.O.C.P - calculate percentage of conserved proteins
##
## A prokaryotic genus can be defined as a group of species
## with all pairwise POCP values higher than 50%
##
## http://jb.asm.org/content/196/12/2210.full

## POCP.The conserved proteins between a pair of genomes were determined by aligning all the protein
# sequences of one genome (query genome) with all the protein sequences of another genome using the
# BLASTP program (20). Proteins from the query genome were considered conserved when they had a BLAST
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

#species1_faa = '/home/hoelzer/RubymineProjects/pocp/Cga_08-1274-3.faa'
#species2_faa = '/home/hoelzer/RubymineProjects/pocp/Cav_10DC88.faa'
#out_dir = '/home/hoelzer/RubymineProjects/pocp/example'
#threads = 6
species1_faa = ARGV[0]
species2_faa = ARGV[1]
out_dir = ARGV[2]
threads = ARGV[3]
#./pocp.rb /mnt/prostlocal2/projects/chlamydia_comparative_study/lisa/prokka/Cav/Cav_10DC88.faa /mnt/prostlocal2/projects/chlamydia_comparative_study/lisa/prokka/Cps/Cps_02DC15.faa ~/pocp_example 6
Pocp.new(species1_faa, species2_faa, out_dir, threads)


