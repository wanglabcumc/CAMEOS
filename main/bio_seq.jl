"""
bio_seq : logic associated with dealing with proteins/translation and
other sequence functions unique to their biological interpretation.
In general our order of aa's is ARNDCQEGHILKMFPSTWYV- .
This is the "alphabetical" order and comes from cross-comptability with Gremlin-deived MRFs.
"""

module bio_seq

using Unicode

islower(s) = all(c->islowercase(c) | isspace(c), s);

nucs = "ACGT"
quick_rev = Dict{String, String}("A" => "T", "C" => "G", "G" => "C", "T" => "A")
codon_table = Dict{String, Char}("..." => '-', "TGT"=>'C',"GAC"=>'D',"TTC"=>'F',"TAG"=>'*',"GTG"=>'V',"CCT"=>'P', "GCT"=>'A',"GGC"=>'G',"CGG"=>'R',
									 "ATT"=>'I',"GAT"=>'D',"CAG"=>'Q',"ATG"=>'M',"CTC"=>'L',"TCT"=>'S',"CGT"=>'R',"ACG"=>'T',"AGA"=>'R',
									 "TGG"=>'W',"TCG"=>'S',"TTA"=>'L',"AGT"=>'S',"CGA"=>'R',"TGC"=>'C',"CAA"=>'Q',"TTG"=>'L',"AAT"=>'N',
									 "AAC"=>'N',"TAA"=>'*',"TAC"=>'Y',"CCA"=>'P',"ACT"=>'T',"TAT"=>'Y',"CAC"=>'H',"CCG"=>'P',"GCA"=>'A',
									 "GAA"=>'E',"ACA"=>'T',"GTA"=>'V',"AGG"=>'R',"AGC"=>'S',"TCA"=>'S',"GGT"=>'G',"GCC"=>'A',"TGA"=>'*',
									 "GTC"=>'V',"TTT"=>'F',"CAT"=>'H',"AAG"=>'K',"AAA"=>'K',"CTT"=>'L',"CTA"=>'L',"ATC"=>'I',"ACC"=>'T',
									 "TCC"=>'S',"CCC"=>'P',"GAG"=>'E',"GCG"=>'A',"CTG"=>'L',"CGC"=>'R',"GTT"=>'V',"GGA"=>'G',"ATA"=>'I',
									 "GGG"=>'G')
aa_to_num = Dict{Char, Int64}('A'=>1, 'R'=>2, 'N'=>3, 'D'=>4, 'C'=>5, 'Q'=>6,
									'E'=>7, 'G'=>8, 'H'=>9, 'I'=>10, 'L'=>11, 'K'=>12, 'M'=>13, 'F'=>14, 'P'=>15,
									'S'=>16, 'T'=>17, 'W'=>18, 'Y'=>19, 'V'=>20, '-'=>21, '*'=>21)

function load_fasta(file_name)
	in_file = open(file_name)
	in_read = readlines(in_file)
	close(in_file)

	seqs = Dict{String, String}()
	cur_head = ""
	cur_seq = String[]
	for line in in_read
		if !(isempty(line)) && line[1] == '>'
			if cur_head != ""
				seqs[cur_head] = join(cur_seq)
			end
			cur_head = strip(line[2:end])
			empty!(cur_seq)
		else
			push!(cur_seq, strip(line))
		end
	end
	seqs[cur_head] = join(cur_seq)
	return seqs
end

function convert_protein(prot_seq)
	prot_len = length(prot_seq)
	big_vector = zeros(Float64, prot_len * 21)
	#and fill it up like a very sparse row vector. One-hot.
	#Order is ARNDCQEGHILKMFPSTWYV-
	prot_fill = Dict{Char, Int64}('A'=>1, 'R'=>2, 'N'=>3, 'D'=>4, 'C'=>5, 'Q'=>6, 'E'=>7, 'G'=>8, 'H'=>9, 'I'=>10,
									'L'=>11, 'K'=>12, 'M'=>13, 'F'=>14, 'P'=>15, 'S'=>16, 'T'=>17, 'W'=>18, 'Y'=>19, 'V'=>20, '-'=>21, 'X'=>21, '*'=>21)
	chr_count = 0
	for chr in prot_seq
		big_vector[(21 * chr_count) + prot_fill[chr]] = 1.0
		chr_count += 1
	end
	return big_vector' #more convenient for it to be a row vector later.
end

function translate(seq, skip=1)
	#adding ellipses for dashes we can't change.
	seq = replace(seq, "*"=>"") #just skip the insertions which we won't change here.
	codons = [seq[x:min(x+2, length(seq))] for x in skip:3:length(seq)]
	out_seq = Char[]
	for codon in codons
		if length(codon) < 3
			continue
		end
		if islower(codon)
			push!(out_seq, lowercase(codon_table[uppercase(codon)]))
		else
			push!(out_seq, codon_table[uppercase(codon)])
		end
	end
	return join(out_seq)
end

function rev_comp(seq)
	return(join([quick_rev[string(B)] for B in reverse(seq)], ""))
end

function get_two_base()
	two_base = [join([one, two]) for one in nucs for two in nucs]
end

function get_three_base() #64.
	three_base = [join([one, two, three]) for one in nucs for two in nucs for three in nucs]
	return three_base
end

function get_four_base() #256 combos.
	four_base = [join([one, two, three, four]) for one in nucs for two in nucs for three in nucs for four in nucs]
	return four_base
end

function get_six_base()
	six_base = [join([one, two, three, four, five, six]) for one in nucs for two in nucs for three in nucs for four in nucs for five in nucs for six in nucs]
	return six_base
end

function precompute_dicodons() #64 * 64 codons = 4096 options in dictionary.
	#Let's precompute all the possible codon pairs to translate very quickly.
	all_bases = get_six_base()
	di_codon_dict = Dict{AbstractString, Tuple{Int64, Int64}}()
	for di_codon in all_bases
		di_codon_dict[di_codon] = translate_to_tuple(di_codon)
	end
	return di_codon_dict
end

function get_base_dict()
	mut_dict = Dict{Int64, Array{AbstractString, 1}}()
	mut_dict[2] = get_two_base()
	mut_dict[3] = get_three_base()
	mut_dict[4] = get_four_base()
	return mut_dict
end

function translate_to_tuple(seq)
	#adding ellipses for dashes we can't change.
	#seq = replace(seq, "*"=>"") #just skip the insertions which we won't change here.
	codons = [seq[x:min(x+2, length(seq))] for x in 1:3:length(seq)]
	out_seq = Char[]
	for codon in codons
		codon = uppercase(codon)
		if length(codon) < 3
			continue
		end
		push!(out_seq, codon_table[uppercase(codon)])
	end
	out_tuple::Tuple{Int64, Int64} = (aa_to_num[out_seq[1]], aa_to_num[out_seq[2]])
	return out_tuple
end

function translate_constrained_maybe_map(nuc, trns; do_map = false)
	gap_poses = Dict{Int64, Int64}()
	for (i, j) in nuc.gap_pos
		if i in keys(gap_poses)
			gap_poses[i] += j
		else
			gap_poses[i] = j
		end
	end

	true_nuc = Char[]
	i_count = 1
	explicit_mapping = Dict{Int64, Int64}()
	aa_counter = 0
	adder = 0

	for unmapped_n in 1:length(nuc.nucs)
		n = unmapped_n + (trns - 1)
		aux_n = n + adder
		if !(n in nuc.skip_trn)
			#then we ignore n for our purposes.
			explicit_mapping[aux_n] = div(aa_counter, 3) + 1
			aa_counter += 1
			push!(true_nuc, nuc.nucs[unmapped_n])
		end
		if n in keys(gap_poses)
			for ii in 1:gap_poses[n]
				for iii in 1:3
					aa_counter += 1
					push!(true_nuc, '.')
				end
			end
		end
	end

	ret_seq = translate(join(true_nuc), 1)
	gap_counter = 0
	if do_map
		return ret_seq, explicit_mapping
	else
		return ret_seq
	end
end

end
