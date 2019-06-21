"""
constants : small set of variables related to translation
"""

aa_table = Dict{String, Array{String, 1}}("A" => ["GCA", "GCC", "GCG", "GCT"], "C" => ["TGT", "TGC"],
            "E" => ["GAG", "GAA"], "D" => ["GAT", "GAC"], "G" => ["GGT", "GGG", "GGA", "GGC"],
            "F" => ["TTT", "TTC"], "I" => ["ATC", "ATA", "ATT"], "H" => ["CAT", "CAC"],
            "K" => ["AAG", "AAA"], "*" => ["TAG", "TGA", "TAA"], "M" => ["ATG"],
            "L" => ["CTT", "CTG", "CTA", "CTC", "TTA", "TTG"], "N" => ["AAC", "AAT"],
            "Q" => ["CAA", "CAG"], "P" => ["CCT", "CCG", "CCA", "CCC"],
            "S" => ["AGC", "AGT", "TCT", "TCG", "TCC", "TCA"],
            "R" => ["AGG", "AGA", "CGA", "CGG", "CGT", "CGC"],
            "T" => ["ACA", "ACG", "ACT", "ACC"], "W" => ["TGG"],
            "V" => ["GTA", "GTC", "GTG", "GTT"], "Y" => ["TAT", "TAC"])

first_letter_codons = Dict{String, Array{String, 1}}("T" => ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TGT", "TGC", "TGG"],
                       "G" => ["GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"],
                       "A" => ["ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG"],
                       "C" => ["CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG"])

last_letter_codons = Dict{String, Array{String, 1}}("C" => ["TTC", "TCC", "TAC", "TGC", "CTC", "CCC", "CAC", "CGC", "ATC", "ACC", "AAC", "AGC", "GTC", "GCC", "GAC", "GGC"],
                      "T" => ["TTT", "TCT", "TAT", "TGT", "CTT", "CCT", "CAT", "CGT", "ATT", "ACT", "AAT", "AGT", "GTT", "GCT", "GAT", "GGT"],
                      "G" => ["TTG", "TCG", "TGG", "CTG", "CCG", "CAG", "CGG", "ATG", "ACG", "AAG", "AGG", "GTG", "GCG", "GAG", "GGG"],
                      "A" => ["TTA", "TCA", "CTA", "CCA", "CAA", "CGA", "ATA", "ACA", "AAA", "AGA", "GTA", "GCA", "GAA", "GGA"])

first_letter_possibilities = Dict{String, Array{String, 1}}("A" => ["I", "M", "T", "N", "K", "S", "R"],
                              "C" => ["L", "P", "H", "Q", "R"],
                              "T" => ["F", "L", "S", "Y", "C", "W"],
                              "G" => ["V", "A", "D", "E", "G"])

top_codon = Dict{Char, AbstractString}('A'=>"GCG", 'C'=>"TGC", 'D'=>"GAC", 'E'=>"GAA", 'F'=>"TTT",
											'G'=>"GGT", 'H'=>"CAT", 'I'=>"ATT", 'K'=>"AAA", 'L'=>"CTG",
											'M'=>"ATG", 'N'=>"AAC", 'P'=>"CCG", 'Q'=>"CAG", 'R'=>"CGT",
											'S'=>"AGC", 'T'=>"ACC", 'V'=>"GTG", 'W'=>"TGG", 'Y'=>"TAT")
