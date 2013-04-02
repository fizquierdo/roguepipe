#!/usr/bin/env ruby

usage = "#{$0} consel.log num_runs"
raise usage unless ARGV.size == 2

# Paese the rank column from the consel log
consellog = ARGV[0]
num_runs = ARGV[1].to_i
raise "consel log #{consellog} not found"  unless File.exist?(consellog)
itemline = false
num_items = num_runs.to_i * 2
ranking = []
File.open(consellog).each_line do |l|
  itemline = true if  l =~ /^# rank item/
  if itemline and  l =~ /^#/ 
    unless l =~ /^# rank item/
      entries = l.split
      ranking << entries[2]
    end
  end
end
puts "Parsed ranking " + ranking.join(',')
puts "Ranking first half " + ranking.slice(0,num_runs).join(',')
puts "Ranking second half " + ranking.slice(num_runs, num_runs).join(',')
# See if the 2 groups have mixed in the ranking
expected_set = (1..num_runs).to_a
best_ranking_set = ranking.slice(0,num_runs)

different_topologies = true
best_ranking_set.each do |item|
  if expected_set.include?(item.to_i)
    #puts "#{item} is the expected set #{expected_set.to_s}"
  else
    different_topologies = false
    #puts "#{item} is not in the expected set #{expected_set.to_s}"
  end
end
# Output a key word to decide if the next iteration is required
if different_topologies
  puts "DIFFERENT_TOPOLOGIES"
else
  puts "SIMILAR_TOPOLOGIES"
end






