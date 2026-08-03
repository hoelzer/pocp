pocp_values = []
pocp_sum = 0
#File.open('pocp-diff.tsv','r').each do |line|
File.open('pocp-diff-old-vs-new-blast.tsv','r').each do |line|
#File.open('brucella-pocp-diff.tsv','r').each do |line|
#File.open('enterococcus-pocp-diff.tsv','r').each do |line|
#File.open('chlamydia-pocp-diff.tsv','r').each do |line|
#File.open('klebsiella-pocp-diff.tsv','r').each do |line|
s = line.split("\t")
    next if line.start_with?('ID')
    line.split("\t").each do |value|
        pocp = value.to_f.abs
        if pocp > 0
            pocp_values.push(pocp)
            pocp_sum += pocp
        end
    end
end
puts "Average percentage difference is #{(pocp_sum / pocp_values.size).round(3)}"
