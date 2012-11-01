require 'rubygems'
require 'fastercsv'
require 'fasta_parser'

outfile = 'results/gene2common_names.csv'
proteome_db_fasta_file = 'data/MOUSE.fasta'

protein_to_common_names = {}
FasterCSV.foreach("data/protein_names.csv") do |row|
	prot_name = row[0]
	protein_to_common_names[prot_name] = nil
end

fasta_dictionary = {}
proteome_db_fap = FastaParser.open(proteome_db_fasta_file)
proteome_db_fap.each do |fasta_entry|
        key1 = fasta_entry.accno
        fasta_entry.desc =~ /(.+)_MOUSE.+GN=(.+)\s/
        key2 = $1
        key3 = $2

        fasta_dictionary[key1] = fasta_entry.desc
        fasta_dictionary[key2] = fasta_entry.desc
        fasta_dictionary[key3] = fasta_entry.desc
end

protein_to_common_names.each_key do |prot_name|
      protein_to_common_names[prot_name] = fasta_dictionary[prot_name]
end

FasterCSV.open(outfile,'w') do |csv|
	csv << %w{ PROTEIN_NAME COMMON_NAME }
        protein_to_common_names.each_key do |prot_name|
              csv << [prot_name, protein_to_common_names[prot_name]]
        end
end
