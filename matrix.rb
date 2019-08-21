## produce a matrix/excel format from the POCP script output

pocps = {}
#strains = ['14-2711_R47', '15-2067_O50', '15-2067_O99'] # the flamingos
strains = ARGV[1].split(",")

Dir.glob("#{ARGV[0]}/*-vs-*").each do |comp|
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

out = File.open("#{ARGV[0]}/#{File.basename(ARGV[0])}.csv",'w')
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
		  pocp = `grep POCP #{ARGV[0]}/#{comp}/pocp.txt | awk 'BEGIN{FS=" "};{print $5}'`.chomp.strip
			out_line += ",#{pocp}"
    end    
  end
  out << out_line << "\n"  
end
out.close

