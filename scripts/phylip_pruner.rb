#!/usr/bin/env ruby

class Phylip
  attr_reader :numtaxa, :seqlen, :seqs
  def initialize(phylipfile)
    raise "File #{phylipfile} does not exist" unless File.exists?(phylipfile)
    @filename = phylipfile
    @seqs = File.open(phylipfile).readlines
    @numtaxa, @seqlen = @seqs[0].split.map{|w| w.to_i}
    @seqs.delete_at(0)
    @seqs.delete_if{|l| l=~ /^\s+$/}
    raise "wrong number of seqs,parsed #{@seqs.size} expected ntaxa #{@numtaxa}" unless @seqs.size == @numtaxa
  end
  def remove_taxa(taxa, pruned_phylip)
    puts "Original size #{@seqs.size}, after removal expect #{@seqs.size - taxa.size}"
    raise "empty list of taxa to prune" if not taxa or taxa.empty?
    taxa.each do |taxon|
      @seqs.delete_if{|l| l.split.first.strip == taxon}
    end
    self.save_as(pruned_phylip)
    puts "Final size #{@seqs.size} saved in #{pruned_phylip}"
  end
  def save_as(newfile)
    self.save_seqs_as(@seqs, newfile)
  end
   def save_seqs_as(seqs, newfile)
    File.open(newfile, "w") do |f|
      f.puts "#{seqs.size} #{@seqlen}"
      seqs.each{|seq| f.puts seq}
    end
  end
end

raise "#{$0} blacklist phylipfile new_phylipfile" unless ARGV.size == 3

blacklist = File.open(ARGV[0]).readlines.to_a.map{|l| l.chomp!}
phylipfile = ARGV[1]
new_phylipfile = ARGV[2]

phy = Phylip.new phylipfile
phy.remove_taxa blacklist, new_phylipfile
