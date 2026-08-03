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

# read in BLAST POCP values
pocp_new = {}
#File.open('blastp/pocp-matrix.tsv','r').each do |line|
File.open('diamond/pocp-matrix.tsv','r').each do |line|
    s = line.split("\t")
    next if line.start_with?('ID')
    species = s[0]
    pocp_new[species] = []
    s.drop(1).each do |pocp|
        pocp_new[species].push(pocp.to_f) 
    end
end

# snaity, they have diff order! need to sort them before diff!
puts pocp_old.keys.join("\t")
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,ctr,cru,pac,sne,wch
puts pocp_new.keys.join("\t")
#cab,cav,cca,cfe,cga,cib,cmu,cpe,cpn,cps,cru,ctr,pac,sne,wch

pocp_sorted = {}
# switch in the arrays entries 10 and 11
pocp_new.each do |species, pocp_values|
    entry_a = pocp_values[10]
    entry_b = pocp_values[11]
    pocp_sorted[species] = pocp_values
    pocp_sorted[species][10] = entry_b
    pocp_sorted[species][11] = entry_a
end
puts pocp_sorted.keys.join("\t")

pocp_combined = {}
species_list.each do |species|
    pocp_values = pocp_sorted[species]
    array = pocp_old[species]
    pos = 0
    pocp_values.each do |pocp|
        if array[pos] == 0
            array[pos] = pocp_values[pos].to_f.round(1)
        end
        pos += 1
    end
    pocp_combined[species] = array
end
puts pocp_combined

# output a new matrix
#output = File.open('blastp/pocp-combined.tsv','w')
output = File.open('diamond/pocp-combined.tsv','w')
output << "ID\t#{species_list.join("\t")}\n"
pocp_combined.each do |species, pocp_array|
    output << "#{species}\t#{pocp_array.join("\t")}\n"
end
output.close



