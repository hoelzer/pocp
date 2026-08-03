# read in old POCP values from 2016 paper
pocp_old = {}
species_list = []
File.open('../old-pocp.csv','r').each do |line|
    s = line.split(",")
    next if line.start_with?('ID')
    species = s[0]
    species_list.push(species)
    pocp_old[species] = []
    s.drop(1).each do |pocp|
        pocp_old[species].push(pocp.to_f) 
    end
end
#puts pocp_old

# read in BLAST POCP values
pocp_blast = {}
File.open('blastp/pocp-matrix.tsv','r').each do |line|
    s = line.split("\t")
    next if line.start_with?('ID')
    species = s[0]
    pocp_blast[species] = []
    s.drop(1).each do |pocp|
        pocp_blast[species].push(pocp.to_f) 
    end
end
puts pocp_blast

# snaity, they have diff order! need to sort them before diff!
puts pocp_old.keys.join(",")
puts pocp_blast.keys.join(",")
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,ctr,cru,pac,sne,wch
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,cru,ctr,pac,sne,wch
pocp_blast_sorted = {}
# switch in the arrays entries 10 and 11
pocp_blast.each do |species, pocp_values|
    entry_a = pocp_values[10]
    entry_b = pocp_values[11]
    pocp_blast_sorted[species] = pocp_values
    pocp_blast_sorted[species][10] = entry_b
    pocp_blast_sorted[species][11] = entry_a
end
puts pocp_blast_sorted

# diff
pocp_diff = {}
species_list.each do |species|
    blast = pocp_blast[species]
    old = pocp_old[species]
    pocp_diff[species] = []
    blast.size.times do |i|
        if old[i] == 0
            pocp_diff[species].push('na')
        else
            pocp_diff[species].push((old[i]-blast[i]).round(2))
        end
    end
end
#puts pocp_diff

# output a new matrix, just output the uper triangle
output = File.open('pocp-diff-old-vs-new-blast.tsv','w')
output << "ID\t#{species_list.join("\t")}\n"
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



