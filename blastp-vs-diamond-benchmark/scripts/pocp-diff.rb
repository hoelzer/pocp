# read in BLAST POCP values
pocp_blast = {}
species_list = []
#File.open('blastp/pocp-matrix.tsv','r').each do |line|
#File.open('results-brucella/blastp/pocp-matrix.tsv','r').each do |line|
#File.open('results-enterococcus/blastp/pocp-matrix.tsv','r').each do |line|
File.open('results-chlamydia/blastp/pocp-matrix.tsv','r').each do |line|
#File.open('results-klebsiella/blastp/pocp-matrix.tsv','r').each do |line|
    s = line.split("\t")
    next if line.start_with?('ID')
    species = s[0]
    species_list.push(species)
    pocp_blast[species] = []
    s.drop(1).each do |pocp|
        pocp_blast[species].push(pocp.to_f) 
    end
end
puts pocp_blast

pocp_diamond = {}
#File.open('diamond/pocp-matrix.tsv','r').each do |line|
#File.open('results-brucella/diamond/pocp-matrix.tsv','r').each do |line|
#File.open('results-enterococcus/diamond/pocp-matrix.tsv','r').each do |line|
File.open('results-chlamydia/diamond/pocp-matrix.tsv','r').each do |line|
#File.open('results-klebsiella/diamond/pocp-matrix.tsv','r').each do |line|
    s = line.split("\t")
    next if line.start_with?('ID')
    species = s[0]
    pocp_diamond[species] = []
    s.drop(1).each do |pocp|
        pocp_diamond[species].push(pocp.to_f) 
    end
end
puts pocp_diamond

# sanity
puts pocp_diamond.keys.join(",")
puts pocp_blast.keys.join(",")
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,cru,ctr,pac,sne,wch
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,cru,ctr,pac,sne,wch
# 
# swap cru and ctr

# switch in the arrays entries 10 and 11
pocp_blast_sorted = {}
pocp_blast.each do |species, pocp_values|
    entry_a = pocp_values[10]
    entry_b = pocp_values[11]
    pocp_blast_sorted[species] = pocp_values
    pocp_blast_sorted[species][10] = entry_b
    pocp_blast_sorted[species][11] = entry_a
end
#puts pocp_blast_sorted
pocp_diamond_sorted = {}
pocp_diamond.each do |species, pocp_values|
    entry_a = pocp_values[10]
    entry_b = pocp_values[11]
    pocp_diamond_sorted[species] = pocp_values
    pocp_diamond_sorted[species][10] = entry_b
    pocp_diamond_sorted[species][11] = entry_a
end
#puts pocp_diamond_sorted

# switch in species list!
species_list_sorted = species_list
entry_a = species_list[10]
entry_b = species_list[11]
species_list_sorted[10] = entry_b
species_list_sorted[11] = entry_a

# diff
pocp_diff = {}
species_list_sorted.each do |species|
    blast = pocp_blast_sorted[species]
    diamond = pocp_diamond_sorted[species]
    pocp_diff[species] = []
    blast.size.times do |i|
        pocp_diff[species].push((blast[i]-diamond[i]).round(3))
    end
end
#puts pocp_diff

# output a new matrix, just output the uper triangle
#output = File.open('pocp-diff.tsv','w')
#output = File.open('brucella-pocp-diff.tsv','w')
#output = File.open('enterococcus-pocp-diff.tsv','w')
output = File.open('chlamydia-pocp-diff.tsv','w')
#output = File.open('klebsiella-pocp-diff.tsv','w')
output << "ID\t#{species_list_sorted.join("\t")}\n"
i = 0
pocp_diff.each do |species, pocp_array|
    pocp_array.shift(i)
    output << "#{species}\t"
    i.times do
        output << "\t"
    end
    output << "#{pocp_array.join("\t")}\n"
    i += 1
end

