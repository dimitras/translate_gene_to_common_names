require 'rubygems'
require 'fastercsv'
require 'fasta_parser'

MW = {
  'A' => 71.03712 ,
  'R' => 156.10112,
  'N' => 114.04293,
  'D' => 115.02695,
  'C' => 103.00919,
  'E' => 129.0426 ,
  'Q' => 128.05858,
  'G' => 57.02147 ,
  'H' => 137.05891,
  'I' => 113.08407,
  'L' => 113.08407,
  'K' => 128.09497,
  'M' => 131.04049,
  'F' => 147.06842,
  'P' => 97.05277 ,
  'S' => 87.03203 ,
  'T' => 101.04768,
  'U' => 150.95364,
  'W' => 186.07932,
  'Y' => 163.06333,
  'V' => 99.06842
}

AVGMASS = {
  'A' => 71.08 ,
  'R' => 156.19,
  'N' => 114.1 ,
  'D' => 115.09,
  'C' => 103.14,
  'E' => 129.12,
  'Q' => 128.13,
  'G' => 57.05 ,
  'H' => 137.14,
  'I' => 113.16,
  'L' => 113.16,
  'K' => 128.17,
  'M' => 131.19,
  'F' => 147.18,
  'P' => 97.12 ,
  'S' => 87.08 ,
  'T' => 101.1 ,
  'U' => 150.03,
  'W' => 186.21,
  'Y' => 163.18,
  'V' => 99.13
}

H2O = 18.015
H = 1.00794


def calc_mw(seq=[])
  mw = 0.0
  seq.each do |aa|
    mw += AVGMASS[aa].to_f
  end
  return mw + H2O
end

# outfile = 'results/proteins_translation.csv'
outfile = 'results/Mouse_AAgenes_uncurated_translation.csv'
all_outfile = 'results/Mouse_AAgenes_uncurated_translation_all.csv'
proteome_db_fasta_file = 'data/MOUSE2013.fasta'
# names_list_file = 'data/protein_names.csv'
names_list_file = 'data/Mouse_AAgenes_uncurated_updated.csv'

genes_to_common_names = Hash.new { |h,k| h[k] = [] }
FasterCSV.foreach(names_list_file) do |row|
  genename = row[0]
  genes_to_common_names[genename] = nil
end

fasta_dictionary = Hash.new { |h,k| h[k] = [] }
proteome_db_fap = FastaParser.open(proteome_db_fasta_file)
proteome_db_fap.each do |fasta_entry|
  isoform = false
  prot_accno = fasta_entry.accno # protein accno
  if prot_accno.include?('-')
    fasta_entry.desc =~ /(.+)_MOUSE.+GN=(.+).*/
    isoform = true
  else
    fasta_entry.desc =~ /(.+)_MOUSE.+GN=(.+)\sPE.+/
  end
  genename = $2 # genename
  mw = calc_mw(fasta_entry.seq.split(''))
  fasta_dictionary[genename] << [prot_accno, fasta_entry.desc, mw, fasta_entry.tag, isoform]
end

genes_to_common_names.each_key do |genename|
  genes_to_common_names[genename] = fasta_dictionary[genename]
end

FasterCSV.open(all_outfile,'w') do |csv|
  csv << %w{ GENENAME PROTEIN_NAME COMMON_NAME MW TAG}
  genes_to_common_names.each do |genename, entries|
    for i in 0..entries.length-1
      csv << [genename, genes_to_common_names[genename][i][0], genes_to_common_names[genename][i][1], genes_to_common_names[genename][i][2], genes_to_common_names[genename][i][3]]
    end
  end
end

# export file with only the non-isoforms and the swiss-prot entries
FasterCSV.open(outfile,'w') do |csv|
  csv << %w{ GENENAME PROTEIN_NAME COMMON_NAME MW}
  genes_to_common_names.each do |genename, entries|
    for i in 0..entries.length-1
      if (genes_to_common_names[genename][i][3] == 'sp') && (genes_to_common_names[genename][i][4] == false)
        csv << [genename, genes_to_common_names[genename][i][0], genes_to_common_names[genename][i][1], genes_to_common_names[genename][i][2]]
      end
    end
  end
end
