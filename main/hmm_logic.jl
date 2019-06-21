"""
hmm_logic : logic associated with dealing with Hidden Markov Models and
associated files, as well as construction of internal versions of this model.
"""

using BioAlignments, BioSymbols, Unicode

islower(s) = all(c->islowercase(c) | isspace(c), s);
#isupper(s) = all(c->isuppercase(c) | isspace(c), s);

align_problem = GlobalAlignment()
align_blosum62 = BLOSUM62;
#Custom parameters that seem to handle "edges" of alignment well: we want a pretty heavy penalty on gaps.
align_model = AffineGapScoreModel(align_blosum62, match=2, mismatch=-1, gap_open=-20, gap_extend=-5)


global INF_NEG
INF_NEG = -1e12

global bases
bases = ["A", "C", "G", "T"]

function get_hmm_len(hmm_file)
	out_stats = read(pipeline(`hmmstat $hmm_file`), String)
	final_line = split(out_stats, "\n")[end-1]
	the_len = split(final_line)[6]
	return parse(Int64, the_len)
end

function state_probs(hmm::types.wholeHMM, node_num::Int)
	if node_num <= size(hmm.state_probs)[1]
		return hmm.state_probs[node_num, 1:end] #this should be 7 things long...
	else
		#Parameters appropriate for end of sequence.
		return [-0.02312, -3.7785, -1e5, -0.61958, -0.77255, 0.0000, -1e5]
	end
end

function emit_match_probs(hmm::types.wholeHMM, node_num::Int, aas::Array{Char, 1})
	return Dict{Char, Float64}(aas[i] => hmm.match_probs[node_num, i] for i in 1:20) #depends on aa, alphabetical. Think it's fine.
end

function emit_insert_probs(hmm::types.wholeHMM, node_num::Int, aas::Array{Char, 1})
	return Dict(aas[i] => hmm.insert_probs[node_num, i] for i in 1:20)
end

function hmm_get_consensus(hmm, hmm_length::Int)
	#This is a simple "return most likely residue at each position" consensus.
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	consensus = ""
	for pos in 1:(hmm_length-2)
		consensus *= string(aas[findmax(hmm.match_probs[pos, 1:end])[2]])
	end
	consensus *= "*"
	return consensus
end

function create_HMM(file_name)
	#Extracts parameters from hmm files. Returns them to construct HMM objects.
	in_file = open(file_name)
	in_read = readlines(in_file)
	close(in_file)

	line_count = 0
	block_mod = 0
	last_line = ""
	skip_block = false

	state_probs = zeros(Float64, 0, 0)
	insert_aa_probs = zeros(Float64, 0, 0)
	match_aa_probs = zeros(Float64, 0, 0)
	hmm_length = 0
	alphabet = ""
	alphabet_length = 20 #for now, will be set again downstream
	collecting = false
	cur_row = -1 #-1 because we start collecting at a distinctive "m->m" line which exists in every hmm file but is ahead of actual parameters.

	for line in in_read[1:end-1]

		if !(collecting)
			if line[1:4] == "LENG"
				hmm_length = parse(Int64, split(line)[2])
				state_probs = zeros(Float64, hmm_length, 7)
				insert_aa_probs = zeros(Float64, hmm_length, alphabet_length)
				match_aa_probs = zeros(Float64, hmm_length, alphabet_length)
			end
		end
		if (collecting)
			section = block_mod % 3
			elements = split(line)

			if section == 0
				cur_row += 1
				hmm_pos = elements[1]
				if hmm_pos == "COMPO"
					skip_block = true
				else
					skip_block = false
					probs = Float64[]
					for ele in elements[2:alphabet_length+1]
						if ele == "*"
							push!(probs, 0.0)
						else
							push!(probs, -parse(Float64, ele))
						end
					end
					match_aa_probs[cur_row, 1:end] = probs
				end

			elseif (section == 1) && !(skip_block)
				elements = split(line)
				probs = Float64[]
				for ele in elements[1:alphabet_length]
					if ele == "*"
						push!(probs, 0.0)
					else
						push!(probs, -parse(Float64, ele))
					end
				end
				insert_aa_probs[cur_row, 1:end] = probs
			elseif (section == 2) && !(skip_block)
				probs = Float64[]
				for ele in elements
					if ele == "*"
						push!(probs, 0.0)
					else
						push!(probs, -parse(Float64, ele))
					end
				end
				state_probs[cur_row, 1:end] = probs
			end

		elseif (occursin("m->m", line)) && (occursin("m->i", line)) && (occursin("i->i", line))
			collecting = true
			alphabet = split(last_line)[2:end]
			alphabet_length = length(alphabet)
		end

		block_mod += 1

		last_line = strip(line)
		line_count += 1
	end

	return hmm_length, state_probs, insert_aa_probs, match_aa_probs
end

function load_hmm_seq(prot_name, hmm, hmm_len)
	prot_dict = Dict{AbstractString, AbstractString}(prot_name => hmm_get_consensus(hmm, hmm_len))
	nucl_dict = Dict{AbstractString, AbstractString}(prot_name => optimize_codons(prot_dict[prot_name]))
	return prot_dict, nucl_dict
end

top_codon = Dict{Char, AbstractString}('A'=>"GCG", 'C'=>"TGC", 'D'=>"GAC", 'E'=>"GAA", 'F'=>"TTT",
																							'G'=>"GGT", 'H'=>"CAT", 'I'=>"ATT", 'K'=>"AAA", 'L'=>"CTG",
																							'M'=>"ATG", 'N'=>"AAC", 'P'=>"CCG", 'Q'=>"CAG", 'R'=>"CGT",
																							'S'=>"AGC", 'T'=>"ACC", 'V'=>"GTG", 'W'=>"TGG", 'Y'=>"TAT")


function get_hmm_consensus(hmm_file) #NOTE: extra dependency (hmmemit)
	hmm_con = read(`hmmemit -c $hmm_file`, String)
	return join(split(hmm_con, '\n')[2:end], "")
end

function get_tails(hmmalign) #HMM alignments sometimes struggle near end of sequence. This is a helper function to collect ends of sequence.
	#We correct the problems with the HMM alignment by just using a standard pairwise aligner.
	cur_pos = 1
	seq_to_date = Char[]
	while islower(hmmalign[cur_pos]) || hmmalign[cur_pos] == '-'
		push!(seq_to_date, hmmalign[cur_pos])
		cur_pos += 1
	end

	end_count = 7
	while end_count > 0 && cur_pos < length(hmmalign)
		if isuppercase(hmmalign[cur_pos])
			push!(seq_to_date, hmmalign[cur_pos])
			cur_pos += 1
			end_count -= 1
		else
			end_count = 0
		end
	end

	out_str = join(seq_to_date)
	return out_str
end

function get_consensus(hmm_file)
	get_hmm_consensus(hmm_file)
end

function get_seqs(seq_keys, gene_name)
	for k in seq_keys
		if occursin("gene:$gene_name ", k)
			return k
		end
	end
	return "N/A"
end

function align_grem(gremlin_prot, hmm_file)
	#This just aligns the HMM consensus with the MRF model. We do this to map between both sequence models.
	#Often the HMM has more insertion points than the MRF since HMMs are better at handling insertions. MRFs we just prune nodes.
	the_str = ">grem\n$gremlin_prot\n"
	hmm_consensus = get_hmm_consensus(hmm_file)
	hmm_out = read(pipeline(`printf "$the_str"`, `hmmalign --outformat A2M $hmm_file -`), String)
	out_hmm = join(split(hmm_out, '\n')[2:end], "")
	out_hmm = fix_hmmalign(out_hmm, hmm_consensus)
	return out_hmm
end

function align_seq(seq_prot, hmm_file, hmm_consensus) #this function is very similar to last, slightly different input parameters.
	the_str = ">seq\n$seq_prot\n"
	hmm_out = read(pipeline(`printf "$the_str"`, `hmmalign --outformat A2M $hmm_file -`), String)
	out_hmm = join(split(hmm_out, '\n')[2:end], "")
	out_hmm = fix_hmmalign(out_hmm, hmm_consensus)
	return out_hmm
end

function redo_alignment(grem_text, sample_nucs, is_mark_not_deg)
	#More logic around fixing issues with model alignments.
	#Here we're trying to map out insertions relative to MRF/HMM.
	k = 0
	skip = Dict{Int64, Array{Int64, 1}}()
	deletions = Dict{Int64, Array{Int64, 1}}()

	for ak in 1:length(sample_nucs) #making a skip/deletion vector for each member of sample_nucs.
		skip[ak] = Int64[]
		deletions[ak] = Int64[]
	end

	for sn in sample_nucs
		k += 1
		nat_path = sn.path[1:end]
		if is_mark_not_deg
			first_aa = sn.path[end][1]
		else
			first_aa = sn.path[end][2]
		end

		for_path = reverse(nat_path) #this is the forward path.
		if is_mark_not_deg
			short_path = [x[1] for x in for_path]
		else
			short_path = [x[2] for x in for_path]
		end

		hmm_pos = 0 #keeps track of HMM-node position. Doesn't increase on lower-case letters.
		grem_pos = 0 #keeps track of MRF-node position. What we are mapping to.
		path_pos = 0
		last_path_node = 0
		last_normal = 1
		while_flag = false

		for char in 1:length(grem_text)
			if !(islower(grem_text[char]))
				hmm_pos += 1
			end

			if !(grem_text[char] == '-')
				grem_pos += 1
			end

			if hmm_pos < minimum(short_path) || hmm_pos > maximum(short_path)
				#Then we are in unaligned zone, we just need to look at what grem/hmm are up to.
				if grem_text[char] == '-'
					push!(skip[k], grem_pos)
				elseif islower(grem_text[char])
					push!(deletions[k], last_normal)
				else
					last_normal = grem_pos
				end

			else
				if !(while_flag)
					path_pos += 1
				end
				while_flag = false

				if path_pos > length(short_path)
					println("Warning: when redoing alignment, HMM_pos is still in progress path is finished.")
				else
					cur_path_pos = short_path[path_pos]
				end

				#Here we're reacting to something specific to sampling trajectory, we deal with it. This is independent of alignment between HMM/GREM.
				while path_pos < length(short_path) && cur_path_pos == last_path_node
					while_flag = true
					#that means CAMEOS stayed at the same hmm node. An insertion.
					#if there is an insertion here from MRF, presumably it would consume the HMM "insertion"
					if islower(grem_text[char])
						#We don't consider this skippable because there is an insertion from the MRF offsetting the HMM insertion.
						grem_pos += 1
						char += 1
					else
						#If just business as usual, though, then we skip the hmm insertion.
						push!(skip[k], grem_pos)
					end
					#Move orward...
					last_path_node = short_path[path_pos]
					last_path_pos = path_pos
					path_pos += 1
					cur_path_pos = short_path[path_pos]
				end

				if !(while_flag)
				#So this is again the territory where not related to skips in path but between HMM/MRF.
					if grem_text[char] == '-'
						#Then an aa in HMM isn't seen in MRF. We skip this.
						#There's a chance it's deleted. In practice this almost never happens.
						push!(skip[k], grem_pos)
					elseif islower(grem_text[char])
						#Then there is an insertion from MRF vs the HMM.
						#Unless there is an insertion from HMM (which would be covered above) this means we have a deletion.
						push!(deletions[k], last_normal)
						path_pos -= 1 #this is not related to path, we need to go back one.
					elseif (short_path[path_pos] - last_path_node > 1) && last_path_node != 0 #because then it's just the first iteration and we don't have a last_path_node.
						for s in 1:(short_path[path_pos] - last_path_node)
							push!(deletions[k], last_normal)
						end
					else
						last_normal = grem_pos
					end
					cur_path_pos = short_path[path_pos]
					last_path_node = short_path[path_pos]
					last_path_pos = path_pos
				end
			end
		end
	end

	return skip, deletions
end

function to_aa(seq)
	return map(x -> convert(AminoAcid, x), collect(strip(seq)))
end

function fix_hmmalign(hmmalign, hmm_consensus)
	#So hmmalign can be very gappy and insertiony near the ends of sequences... so we want to do a global pairwise alignment on it.
	#We are now going to use Bio.jl's pairalign() function.
	pure_seq = uppercase(replace(hmmalign, "-"=>""))

	start_seq = get_tails(hmmalign)
	original_up_len = length(start_seq)
	edited_start = uppercase(replace(start_seq, "-"=>""))
	first_aligned = findfirst(z -> isuppercase(z), collect(start_seq)) #inserts before an aligned char are not inserts, they are unaligned.

	up_len = length(start_seq) - count(z -> islower(z), start_seq)
	rev_seq = reverse(get_tails(reverse(hmmalign)))
	original_down_len = length(rev_seq)
	edited_end = uppercase(replace(rev_seq, "-"=>""))
	down_len = length(rev_seq) - count(z -> islower(z), rev_seq)

	head_seq = uppercase(string(hmm_consensus[1:up_len]))
	tail_seq = uppercase(string(hmm_consensus[end - down_len + 1: end]))

	#align_problem is a global function that is just defining we want a global alignment here.

	head_align = alignment(pairalign(align_problem, to_aa(edited_start),
							to_aa(head_seq), align_model))
	tail_align = alignment(pairalign(align_problem, to_aa(edited_end),
							to_aa(tail_seq), align_model))

	head_ali_top = collect(join([x for (x, _) in head_align]))
	head_ali_bot = collect(join([x for (_, x) in head_align]))

	for head_count in 1:length(head_ali_top)
		if head_ali_bot[head_count] == '-'
			head_ali_top[head_count] = lowercase(head_ali_top[head_count])
		end
	end
	head_ali_top = join(head_ali_top)

	tail_ali_top = collect(join([x for (x, _) in tail_align]))
	tail_ali_bot = collect(join([x for (_, x) in tail_align]))

	for tail_count in 1:length(tail_ali_top)
		if tail_ali_bot[tail_count] == '-'
			tail_ali_top[tail_count] = lowercase(tail_ali_top[tail_count])
		end
	end
	tail_ali_top = join(tail_ali_top)

	return head_ali_top * hmmalign[original_up_len + 1 : (end - original_down_len)] * tail_ali_top
end

function hmm_name(gene_name)
	hmm_path_file = open("hmkt_gene_table.txt")
	hmm_read = readlines(hmm_path_file)
	close(hmm_path_file)
	hmm_read = filter(x -> !isempty(x), hmm_read)
	for line in hmm_read
		if '\t' in line
			spl = split(strip(line), '\t')
			if spl[1] == gene_name
				return strip(spl[4])
			end
		end
	end
	return ""
end

#Basically we also find it useful to trace the protein sequence to the hmm.
#We want wild-type sequence for non-overlapping regions, and we want to align it to hmm because hmm is sort of our "base positioning"
#of the models we're constructing.
function gen_hmm_trace(gene_name, gene_hmm)
	cds_seqs = bio_seq.load_fasta("cds.fasta")
	prot_seqs = bio_seq.load_fasta("proteins.fasta")

	local cds_seq, prot_seq
	cds_seq = cds_seqs[gene_name]
	prot_seq = prot_seqs[gene_name]

	hmm_consensus = get_consensus(gene_hmm)
	hmm_align = align_seq(prot_seq, gene_hmm, hmm_consensus)

	i_count = 1
	hmm_counter = 1
	to_delete = Int64[]
	to_add = Dict{Int64, Char}()

	for i in hmm_align
		if i == '-'
			to_add[i_count] = hmm_consensus[hmm_counter]
		end
		if islower(i)
			push!(to_delete, i_count)
		else
			hmm_counter += 1
		end
		i_count += 1
	end

	final_aas = Char[]
	final_codons = AbstractString[]

	aa_count = 1
	cds_count = 1
	for aa in hmm_align #hmm_align is a protein sequence, aligned.
		if aa_count in to_delete
			#Then we don't want to include this codon in subsequent step.
			cds_count += 1
		else
			if !(aa_count in keys(to_add)) #this is the normal case.
				push!(final_aas, aa)
				push!(final_codons, cds_seq[((cds_count - 1) * 3) + 1: cds_count * 3])
				cds_count += 1
			elseif aa_count in keys(to_add)
				push!(final_aas, to_add[aa_count])
				push!(final_codons, top_codon[to_add[aa_count]])
				#we don't iterate on cds_count because we haven't added natural codons.
			end
		end
		aa_count += 1 #aa_count has more to do with alignment.
	end

	final_cds = join(final_codons)
	final_prot_seq = join(final_aas)
	return final_cds, final_prot_seq
end

function optimize_codons(sequence::AbstractString)
	#most popular codons taken from first table at http://openwetware.org/wiki/Escherichia_coli/Codon_usage [super simple way to find codons, not claiming it's very smart]
	most_popular = Dict{Char, AbstractString}('G' => "GGC", 'E' => "GAA", 'D' => "GAT", 'V' => "GTG", 'A' => "GCG",
                  'R' => "CGC", 'K' => "AAA", 'N' => "AAC", 'M' => "ATG", 'I' => "ATT",
                  'T' => "ACC", 'W' => "TGG", 'C' => "TGC", '*' => "TAA", 'F' => "TTT",
                  'S' => "AGC", 'Q' => "CAG", 'H' => "CAT", 'L' => "CTG", 'P' => "CCG", 'Y'=>"TAT")
	return join([most_popular[aa] for aa in sequence])
end

function top_n(el_tuple, val_keys)
	a_max = findmax(el_tuple)
	max_key::Int = val_keys[a_max[2]]
	max_value = a_max[1]
	return (max_key, max_value)
end

### All of these give_log() functions look at the HMM probabilities of regions *around* the overlap.
#The idea is that this allows us to calculate full-length HMM probabilities, not just probabilities for overlap region.
function follow_prot_give_log_2d(hmm, prot, limit)
	sum_logs = Float32[0]
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for position in 1:(limit-2)
		push!(sum_logs, sum_logs[position] + state_probs(hmm, position)[1] + emit_match_probs(hmm, position, aas)[prot[position]])
	end
	if limit - 1 > 0
		push!(sum_logs, sum_logs[limit-1] + state_probs(hmm, limit -1)[1])
	end
	return sum_logs * ones(Float32, 1, 64) #I think this makes a big matrix of this vector.
end

function follow_prot_give_log(hmm, prot, limit)
	sum_logs = Float64[0]
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for position in 1:(limit-2)
		push!(sum_logs, sum_logs[position] + state_probs(hmm, position)[1] + emit_match_probs(hmm, position, aas)[prot[position]])
	end
	if limit - 1 > 0
		push!(sum_logs, sum_logs[limit-1] + state_probs(hmm, limit -1)[1])
	end
	return sum_logs
end

function tail_prot_give_log_2d(hmm, prot, limit)
	sum_logs = Float64[0]
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for position in 1:(limit-2)
		push!(sum_logs, sum_logs[position] + state_probs(hmm, position)[1] + emit_match_probs(hmm, position, aas)[uppercase(prot[position])])
	end
	push!(sum_logs, sum_logs[limit-1] + state_probs(hmm, limit - 1)[1])
	total = sum_logs[end]
	return [total - sum_logs[y] for y in 1:length(sum_logs)]
end

function tail_prot_give_log(hmm, prot, limit)
	sum_logs = Float64[0]
	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for position in 1:(limit-2)
		push!(sum_logs, sum_logs[position] + state_probs(hmm, position)[1] + emit_match_probs(hmm, position, aas)[prot[position]])
	end
	push!(sum_logs, sum_logs[limit-1] + state_probs(hmm, limit - 1)[1])
	total = sum_logs[end]
	return [total - sum_logs[y] for y in 1:length(sum_logs)]
end
