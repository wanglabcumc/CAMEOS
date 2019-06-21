"""
dyn_frwd takes care of "forward" construction of double-encoding sequences during the HMM-parameterized step of CAMEOS.
	Orchestrates the construction of the tensor containing appropriate scores for encodings across positions in gene x and y.
	Calculates the costs of substitutions and transitions between states and is the place where recursion back to the
	start of the sequence (through code in dyn_back) is initiated. "dyn" from "Dyn"amic programming.
"""

include("math.jl")
include("dyn_back.jl")

FORCE_MET = true #whether we need all proteins to start with a M. Generally true, unless dealing with fragments.
trinucs = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG",
				 "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC",
				 "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA",
				 "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]

function recurse_back(X_hmm, Y_hmm, short_node, far_node, cur_state, cur_score, last_letter, full_length, cumulative_score, path, bbt, thresh=1.0, temperature=1.0)

	short_mm::Float64, short_md::Float64, short_im::Float64,
	short_dm::Float64, short_dd::Float64 = state_probs(X_hmm, short_node - 1)[[1, 3, 4, 6, 7]]

	far_mm::Float64, far_md::Float64, far_im::Float64,
	far_dm::Float64, far_dd::Float64 = state_probs(Y_hmm, far_node - 1)[[1, 3, 4, 6, 7]]

	short_mi::Float64, short_ii::Float64 = state_probs(X_hmm, short_node)[[2, 5]]
	far_mi::Float64, far_ii::Float64 = state_probs(Y_hmm, far_node)[[2, 5]]

	all_vals::Array{Float64, 1} = [short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii]

	next_state = 0
	next_codon = ""
	next_stop = 0

	# 1) MxMy, 2) IxMy, 3) DxMy, 4) MxIy, 5) DxIy, 6) MxDy, 7) IxDy
	if cur_state == 1 #Meaning we're in MxMy
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.mxmy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 2
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.ixmy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 3
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.dxmy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 4
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.mxiy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 5
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.dxiy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 6
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.mxdy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	elseif cur_state == 7
		next_state, next_letter, next_codon, cur_score, cumulative_score = dyn_back.ixdy(trinucs, last_letter, short_node, far_node, cumulative_score, cur_score, all_vals, bbt, thresh, temperature)
	end

	#Depending on the state, we move zero or one aa closer to start in both frames.
	if next_state == 1 || next_state == 3 || next_state == 4 || next_state == 5 || next_state == 6
		short_node -= 1
	end
	if next_state == 1 || next_state == 2 || next_state == 3 || next_state == 6 || next_state == 7
		far_node -= 1
	end

	push!(path, (short_node, far_node, next_state))
	if short_node == 1 || far_node == 1 #then we're at a base case of recursion being done.
		return full_length, path, cumulative_score, cur_score
	else
		new_full_length = next_codon * full_length
		return recurse_back(X_hmm, Y_hmm, short_node, far_node, next_state, cur_score, string("ACGT"[next_letter]), new_full_length, cumulative_score, path, bbt, thresh, temperature)
	end
end

function state_transition_calcs(tensor::Array{Float32, 4}, short_node::Int, far_node::Int, X_hmm::types.wholeHMM, Y_hmm::types.wholeHMM, the_keys::Array{Int, 1}, the_vals::Array{Float64, 1},
                      mxmy_calcs::Array{Float64, 1}, mxiy_calcs::Array{Float64, 1}, mxdy_calcs::Array{Float64, 1}, ixmy_calcs::Array{Float64, 1},
                      dxmy_calcs::Array{Float64, 1}, dxiy_calcs::Array{Float64, 1}, ixdy_calcs::Array{Float64, 1},
                      mxmy_keys::Array{Int, 1}, mxiy_keys::Array{Int, 1}, mxdy_keys::Array{Int, 1}, ixmy_keys::Array{Int, 1},
                      dxmy_keys::Array{Int, 1}, ixdy_keys::Array{Int, 1}, dxiy_keys::Array{Int, 1})
	mxmy_back = tensor[short_node - 1, far_node - 1, 1:7, 1:64]
	mxiy_back = tensor[short_node - 1, far_node, 1:5, 1:64]
	mxdy_back = tensor[short_node - 1, far_node - 1, [1, 2, 3, 6, 7], 1:64]
	ixmy_back = tensor[short_node, far_node - 1, [1, 2, 4, 6, 7], 1:64]
	dxmy_back = tensor[short_node - 1, far_node - 1, [1, 3, 4, 5, 6], 1:64]
	ixdy_back = tensor[short_node, far_node - 1, [1, 2, 6, 7], 1:64]
	dxiy_back = tensor[short_node - 1, far_node, [1, 3, 4, 5], 1:64]

	short_mm::Float64, short_md::Float64, short_im::Float64,
	short_dm::Float64, short_dd::Float64 = state_probs(X_hmm, short_node - 1)[[1, 3, 4, 6, 7]]

	far_mm::Float64, far_md::Float64, far_im::Float64,
	far_dm::Float64, far_dd::Float64 = state_probs(Y_hmm, far_node - 1)[[1, 3, 4, 6, 7]]

	short_mi::Float64, short_ii::Float64 = state_probs(X_hmm, short_node)[[2, 5]]
	far_mi::Float64, far_ii::Float64 = state_probs(Y_hmm, far_node)[[2, 5]]

	#In all these cases we are finding the best option for the next state, and adding the cost of the state transition as well.
	mxmy_calcs[1] = maximum(mxmy_back[1, 1:end] .+ (short_mm + far_mm))
	mxmy_calcs[2] = maximum(mxmy_back[2, 1:end] .+ (short_im + far_mm))
	mxmy_calcs[3] = maximum(mxmy_back[3, 1:end] .+ (short_dm + far_mm))
	mxmy_calcs[4] = maximum(mxmy_back[4, 1:end] .+ (short_mm + far_im))
	mxmy_calcs[5] = maximum(mxmy_back[5, 1:end] .+ (short_dm + far_im))
	mxmy_calcs[6] = maximum(mxmy_back[6, 1:end] .+ (short_mm + far_dm))
	mxmy_calcs[7] = maximum(mxmy_back[7, 1:end] .+ (short_im + far_dm))
	the_keys[1], the_vals[1] = top_n(mxmy_calcs, mxmy_keys) #these are all correct.

	ixmy_calcs[1] = maximum(ixmy_back[1, 1:end] .+ (short_mi + far_mm))
	ixmy_calcs[2] = maximum(ixmy_back[2, 1:end] .+ (short_ii + far_mm))
	ixmy_calcs[3] = maximum(ixmy_back[3, 1:end] .+ (short_mi + far_im))
	ixmy_calcs[4] = maximum(ixmy_back[4, 1:end] .+ (short_mi + far_dm))
	ixmy_calcs[5] = maximum(ixmy_back[5, 1:end] .+ (short_ii + far_dm))
	the_keys[2], the_vals[2] = top_n(ixmy_calcs, ixmy_keys)

	dxmy_calcs[1]= maximum(dxmy_back[1, 1:end] .+ (short_md + far_mm))
	dxmy_calcs[2] = maximum(dxmy_back[2, 1:end] .+ (short_dd + far_mm))
	dxmy_calcs[3] = maximum(dxmy_back[3, 1:end] .+ (short_md + far_im))
	dxmy_calcs[4] = maximum(dxmy_back[4, 1:end] .+ (short_dd + far_im))
	dxmy_calcs[5] = maximum(dxmy_back[5, 1:end] .+ (short_md + far_dm))
	the_keys[3], the_vals[3] = top_n(dxmy_calcs, dxmy_keys)

	mxiy_calcs[1] = maximum(mxiy_back[1, 1:end] .+ (short_mm + far_mi))
	mxiy_calcs[2] = maximum(mxiy_back[2, 1:end] .+ (short_im + far_mi))
	mxiy_calcs[3] = maximum(mxiy_back[3, 1:end] .+ (short_dm + far_mi))
	mxiy_calcs[4] = maximum(mxiy_back[4, 1:end] .+ (short_mm + far_ii))
	mxiy_calcs[5] = maximum(mxiy_back[5, 1:end] .+ (short_dm + far_ii))
	the_keys[4], the_vals[4] = top_n(mxiy_calcs, mxiy_keys)

	dxiy_calcs[1] = maximum(dxiy_back[1, 1:end] .+ (short_md + far_mi))
	dxiy_calcs[2] = maximum(dxiy_back[2, 1:end] .+ (short_dd + far_mi))
	dxiy_calcs[3] = maximum(dxiy_back[3, 1:end] .+ (short_md + far_ii))
	dxiy_calcs[4] = maximum(dxiy_back[4, 1:end] .+ (short_dd + far_ii))
	the_keys[5], the_vals[5] = top_n(dxiy_calcs, dxiy_keys)

	mxdy_calcs[1] = maximum(mxdy_back[1, 1:end] .+ (short_mm + far_md))
	mxdy_calcs[2] = maximum(mxdy_back[2, 1:end] .+ (short_im + far_md))
	mxdy_calcs[3] = maximum(mxdy_back[3, 1:end] .+ (short_dm + far_md))
	mxdy_calcs[4] = maximum(mxdy_back[4, 1:end] .+ (short_mm + far_dd))
	mxdy_calcs[5] = maximum(mxdy_back[5, 1:end] .+ (short_im + far_dd))
	the_keys[6], the_vals[6] = top_n(mxdy_calcs, mxdy_keys)

	ixdy_calcs[1] = maximum(ixdy_back[1, 1:end] .+ (short_mi + far_md))
	ixdy_calcs[2] = maximum(ixdy_back[2, 1:end] .+ (short_ii + far_md))
	ixdy_calcs[3] = maximum(ixdy_back[3, 1:end] .+ (short_mi + far_dd))
	ixdy_calcs[4] = maximum(ixdy_back[4, 1:end] .+ (short_ii + far_dd))
	the_keys[7], the_vals[7] = top_n(ixdy_calcs, ixdy_keys)

end

function opt_tensor_core(short_hmm::types.wholeHMM, short_node::Int, far_hmm::types.wholeHMM, far_node::Int,
							assumed_letter::Int, next_letter::Int,
							short_aa_match_probs::Dict{Char, Float64}, short_aa_insert_probs::Dict{Char, Float64},
						 	far_aa_match_probs::Dict{Char, Float64}, far_aa_insert_probs::Dict{Char, Float64},
							aa7_score::Array{Float64, 2},
							codon_array::Array{Char, 1}, codon_array_rev::Array{Char, 1},
							short_aa::Char, far_aa::Char, short_aa_match::Float64,
							far_aa_match::Float64, short_aa_insert::Float64, far_aa_insert::Float64, frame::AbstractString, max_i::Int, max_j::Int)

	#for dinuc in dinucs
	flag = false
	for dinuc in 1:16
		short_aa = codon_array[assumed_letter + dinuc]
		far_aa = codon_array[(dinuc - 1)*4 + next_letter]

		if (short_node < max_i + 0)
			if short_aa == '*'
				continue
			end
			short_aa_match = short_aa_match_probs[short_aa]
			short_aa_insert = short_aa_insert_probs[short_aa]
		else
			if short_aa != '*'
				continue
			end
			short_aa_match = -0.01
			short_aa_insert = -1.5
		end

		if (far_node < max_j + 0)
			if far_aa == '*'
				continue
			end
			far_aa_match = far_aa_match_probs[far_aa]
			far_aa_insert = far_aa_insert_probs[far_aa]
		else
			if far_aa != '*'
				continue
			end
			far_aa_match = -0.01
			far_aa_insert = -1.5
		end

		if FORCE_MET && short_node == 2 && short_aa != 'M'
			continue
		end

		if FORCE_MET && far_node == 2 && far_aa != 'M'
			continue
		end #just don't consider these options... should only invoke assumed_letter is 'A' otherwise no option will work here.

		flag = true

		aa7_score[1, dinuc] = short_aa_match + far_aa_match
		aa7_score[2, dinuc] = short_aa_insert + far_aa_match
		aa7_score[3, dinuc] = short_aa_insert + far_aa_match
		aa7_score[4, dinuc] = short_aa_match + far_aa_insert
		aa7_score[5, dinuc] = short_aa_insert + far_aa_insert
		aa7_score[6, dinuc] = short_aa_match + far_aa_insert
		aa7_score[7, dinuc] = short_aa_insert + far_aa_insert
	end
	return flag
end

function initialize_4d(X_len::Int, Y_len::Int, x_init::Array{Float32, 2}, y_init::Array{Float32, 2}, min_over=50)
	tensor::Array{Float32, 4} = zeros(Float32, X_len, Y_len, 7, 64)
	tensor[1, (Y_len - min_over + 1):end, 1, 1:end] = fill(-1e12, 1, min_over, 64)
	tensor[(X_len - min_over + 1):end, 1, 1, 1:end] = fill(-1e12, min_over, 1, 64)
	if min_over != Y_len
		tensor[1, 1:(Y_len - min_over), 1, 1:end] = x_init
	end
	if min_over != X_len
		tensor[1:(X_len - min_over), 1, 1, 1:end] = y_init
	end
	tensor[1, 1:end, 2:7, 1:end] = fill(-1e12, 1, Y_len, 6, 64)
	tensor[1:end, 1, 2:7, 1:end] = fill(-1e12, X_len, 1, 6, 64)
	return tensor
end

function fill_in_tensor(X_hmm::types.wholeHMM, Y_hmm::types.wholeHMM, X_len::Int, Y_len::Int, x_init::Array{Float32, 2}, y_init::Array{Float32, 2}, X_prot::AbstractString, Y_prot::AbstractString, min_len::Int, frame::AbstractString)
	letter_tensors = Dict{AbstractString, Array{Float32, 4}}("A" => initialize_4d(X_len, Y_len, x_init, y_init, min_len),
										 "C" => initialize_4d(X_len, Y_len, x_init, y_init, min_len),
										 "G" => initialize_4d(X_len, Y_len, x_init, y_init, min_len),
										 "T" => initialize_4d(X_len, Y_len, x_init, y_init, min_len))

	maximal_scores = fill(-1.0e12, 7, 64, 4)

	codon_table = Dict{AbstractString, Char}("TGT"=>'C',"GAC"=>'D',"TTC"=>'F',"TAG"=>'*',"GTG"=>'V',"CCT"=>'P', "GCT"=>'A',"GGC"=>'G',"CGG"=>'R',
										 "ATT"=>'I',"GAT"=>'D',"CAG"=>'Q',"ATG"=>'M',"CTC"=>'L',"TCT"=>'S',"CGT"=>'R',"ACG"=>'T',"AGA"=>'R',
										 "TGG"=>'W',"TCG"=>'S',"TTA"=>'L',"AGT"=>'S',"CGA"=>'R',"TGC"=>'C',"CAA"=>'Q',"TTG"=>'L',"AAT"=>'N',
										 "AAC"=>'N',"TAA"=>'*',"TAC"=>'Y',"CCA"=>'P',"ACT"=>'T',"TAT"=>'Y',"CAC"=>'H',"CCG"=>'P',"GCA"=>'A',
										 "GAA"=>'E',"ACA"=>'T',"GTA"=>'V',"AGG"=>'R',"AGC"=>'S',"TCA"=>'S',"GGT"=>'G',"GCC"=>'A',"TGA"=>'*',
										 "GTC"=>'V',"TTT"=>'F',"CAT"=>'H',"AAG"=>'K',"AAA"=>'K',"CTT"=>'L',"CTA"=>'L',"ATC"=>'I',"ACC"=>'T',
										 "TCC"=>'S',"CCC"=>'P',"GAG"=>'E',"GCG"=>'A',"CTG"=>'L',"CGC"=>'R',"GTT"=>'V',"GGA"=>'G',"ATA"=>'I',
										 "GGG"=>'G')

	codon_array::Array{Char, 1} =	['K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P',
										 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G',
										 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F']

	codon_array_rev::Array{Char, 1} = ['F', 'V', 'L', 'I', 'C', 'G', 'R', 'S', 'S', 'A', 'P', 'T', 'Y', 'D', 'H', 'N', 'L', 'V', 'L', 'M', 'W', 'G',
											'R', 'R', 'S', 'A', 'P', 'T', '*', 'E', 'Q', 'K', 'F', 'V', 'L', 'I', 'C', 'G', 'R', 'S', 'S', 'A', 'P', 'T',
											'Y', 'D', 'H', 'N', 'L', 'V', 'L', 'I', '*', 'G', 'R', 'R', 'S', 'A', 'P', 'T', '*', 'E', 'Q', 'K']

	dinucs = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

	mxmy_keys::Array{Int, 1} = [111, 211, 311, 411, 511, 611, 711]
	mxiy_keys::Array{Int, 1} = [112, 212, 312, 412, 512]
	mxdy_keys::Array{Int, 1} = [113, 213, 313, 613, 713]
	ixmy_keys::Array{Int, 1} = [121, 221, 421, 621, 721]
	dxmy_keys::Array{Int, 1} = [131, 331, 431, 531, 631]
	ixdy_keys::Array{Int, 1} = [123, 223, 623, 723]
	dxiy_keys::Array{Int, 1} = [132, 332, 432, 532]

	aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
				 'T', 'V', 'W', 'Y']

	short_aa::Char = 'A'
	far_aa::Char = 'A'
	short_aa_match::Float64 = 0.0
	far_aa_match::Float64 = 0.0
	short_aa_insert::Float64 = 0.0
	far_aa_insert::Float64 = 0.0

	base_counter::Int = 1
	set_up_counter::Int = 0

	mxmy_calcs::Array{Float64, 1} = zeros(7)
	mxiy_calcs::Array{Float64, 1} = zeros(5)
	mxdy_calcs::Array{Float64, 1} = zeros(5)
	ixmy_calcs::Array{Float64, 1} = zeros(5)
	dxmy_calcs::Array{Float64, 1} = zeros(5)
	dxiy_calcs::Array{Float64, 1} = zeros(4)
	ixdy_calcs::Array{Float64, 1} = zeros(4)

	the_keys::Array{Int, 1} = zeros(Int, 7)
	the_vals::Array{Float64, 1} = zeros(Float64, 7)

	aa7_scores::Array{Float64, 2} = fill(-1000.0, 7, 16)

	min_i::Int = 2
	max_i::Int = X_len
	min_j::Int = 2
	max_j::Int = Y_len

	for i in min_i:max_i
		short_node = i

		if short_node < X_len
			short_aa_match_probs::Dict{Char, Float64} = emit_match_probs(X_hmm, short_node - 1, aas)
			short_aa_insert_probs::Dict{Char, Float64} = emit_insert_probs(X_hmm, short_node - 1, aas)
		else
			short_aa_match_probs = Dict{Char, Float64}(aa => -1e10 for aa in aas)
			short_aa_insert_probs = Dict{Char, Float64}(aa => -1e10 for aa in aas)
		end

		for j in min_j:max_j
			did_set = false
			far_node = j

			if far_node < Y_len
				far_aa_match_probs::Dict{Char, Float64} = emit_match_probs(Y_hmm, far_node - 1, aas)
				far_aa_insert_probs::Dict{Char, Float64} = emit_insert_probs(Y_hmm, far_node - 1, aas)
			else
				far_aa_match_probs = Dict{Char, Float64}(aa => -1e10 for aa in aas)
				far_aa_insert_probs = Dict{Char, Float64}(aa => -1e10 for aa in aas)
			end

			maximal_scores = fill(-1.0e12, 7, 64, 4)

			set_up_counter = 0

			for set_up in bases
				#This basically fills up the values of the_keys and the_vals with top scores for transitions to every other state.
				state_transition_calcs(letter_tensors[set_up], i, j, X_hmm, Y_hmm, the_keys, the_vals,
										 mxmy_calcs, mxiy_calcs, mxdy_calcs, ixmy_calcs, dxmy_calcs, dxiy_calcs, ixdy_calcs,
										 mxmy_keys, mxiy_keys, mxdy_keys, ixmy_keys, dxmy_keys, ixdy_keys, dxiy_keys)

				base_counter = 0
				for next_letter in 1:4
					base_counter += 1

					aa7_scores = fill(-100000.0, 7, 16)
					#Scores each of the states.
					fxn_call = opt_tensor_core(X_hmm, short_node, Y_hmm, far_node, 16 * set_up_counter, base_counter,
												short_aa_match_probs, short_aa_insert_probs, far_aa_match_probs, far_aa_insert_probs,
												aa7_scores, codon_array, codon_array_rev, short_aa, far_aa,
												short_aa_match, far_aa_match, short_aa_insert, far_aa_insert, frame, max_i, max_j)

					total_score = broadcast(+, the_vals, aa7_scores) #So this is a 7x16 matrix.

					maximal_scores[1:end, (1:16) .+ (next_letter-1)*16, set_up_counter + 1] = total_score
				end
				set_up_counter += 1
			end
			next_letter = 0
			for mat_letter in bases
				next_letter += 1
				for mat_num in 1:7
					options = vec(maximal_scores[mat_num, (1:16) .+ (next_letter-1)*16, 1:4])
					letter_tensors[mat_letter][i, j, mat_num, 1:end] = options[1:64]
				end
			end
		end
	end
	return letter_tensors
end

function sample_start_point_on_range(X_hmm, X_prot, X_len, Y_hmm, Y_prot, Y_len, num_samples, X_range, Y_range, thresh=1.0, soft_temp=1.0, first_softmax_temp=0.01, gen_bbt=false)
	#Go along all endpoints in all letters and probabilistically choose a start point based on scores.
	#NOW WITH X_RANGE and Y_RANGE TO INDICATE ACCEPTABLE PARTS OF MATRIX FROM WHICH TO SAMPLE.
	vecs_to_sample_x = Dict{AbstractString, Array{Float64, 1}}() #is the score a Float64?
	vecs_to_sample_y = Dict{AbstractString, Array{Float64, 1}}()
	letter_legal_x = Dict{AbstractString, Array{Float64, 1}}()
	letter_legal_y = Dict{AbstractString, Array{Float64, 1}}()
	best_letter_count = 0

	tail_offset = 0.0

	term_x = false
	term_y = false
	if length(X_range) == 1 && length(Y_range) > 1
		term_x = true #this means we want to go from last letter of x backward.
	elseif length(Y_range) == 1 && length(X_range) > 1
		term_y = true
	else
		return false #this will raise an error.
	end

	if !(gen_bbt == false)
		bbt = gen_bbt
	end

	MIN_STOP_OFFSET = 2

	L_values = Dict{Any, Any}()

	for letter_iterator in ["A", "C", "G", "T"]
		top_L = bbt[letter_iterator][1:end, 1:end, 1:end, 1:end]
		best_letter_count += 1
		for integer in 1:64
			top_L[1:X_len, Y_len:Y_len, 1:7, integer:integer] += reshape(repeat(tail_prot_give_log_2d(X_hmm, X_prot, X_len), 1, 7), X_len, 1, 7)
			top_L[X_len:X_len, 1:Y_len, 1:7, integer:integer] += reshape(repeat(tail_prot_give_log_2d(Y_hmm, Y_prot, Y_len), 1, 7), 1, Y_len, 7)
		end
		L_values[letter_iterator] = top_L[1:end, 1:end, 1:end, 1:end]

		final_x_points = top_L[X_len, Y_range, 1:7, 1:64]
		final_y_points = top_L[X_range, Y_len, 1:7, 1:64]

		vec_to_sample_x = vec(final_x_points)
		vec_to_sample_y = vec(final_y_points)
		vecs_to_sample_x[letter_iterator] = vec_to_sample_x
		vecs_to_sample_y[letter_iterator] = vec_to_sample_y
	end

	ADJUSTED_Y_LEN = length(Y_range)
	ADJUSTED_X_LEN = length(X_range)

	all_scores = Float64[]
	for letter_iterator in ["A", "C", "G", "T"]
		for i in vecs_to_sample_x[letter_iterator]
			if term_x
				push!(all_scores, i)
			end
		end
		for j in vecs_to_sample_y[letter_iterator]
			if term_y
				push!(all_scores, j)
			end
		end
	end

	first_steps = math.softmax(all_scores, num_samples, first_softmax_temp)
	if term_x
		sample_x_len = length(vecs_to_sample_x["A"])
		sample_y_len = 0
	elseif term_y
		sample_x_len = 0
		sample_y_len = length(vecs_to_sample_y["A"])
	end
	per_letter = sample_x_len + sample_y_len

	results = Any[]
	first_step_count = 1
	for first_step in first_steps
		first_step_count += 1

		letter_value = div(first_step - 1, per_letter) + 1
		nuc_start = string("ACGT"[letter_value])
		pos_in = ((first_step - 1) % per_letter) + 1
		alt_aa_pos = 0
		is_x = pos_in < sample_x_len
		if is_x
			aa_pos = (Y_range[1] - 1) + (((pos_in - 1) % (ADJUSTED_Y_LEN)) + 1) #maybe -1 on the Y_range business...
			alt_aa_pos = ((pos_in - 1) % (ADJUSTED_Y_LEN)) + 1
			cur_state = div((pos_in - 1) % ((ADJUSTED_Y_LEN) * 7), (ADJUSTED_Y_LEN)) + 1
			last_codon = div(pos_in - 1, ((ADJUSTED_Y_LEN) * 7)) + 1
		else
			transformed = mod((first_step - 1), per_letter) + 1 - sample_x_len
			aa_pos = (X_range[1] - 1) + (((transformed - 1) % (ADJUSTED_X_LEN)) + 1)
			alt_aa_pos = ((transformed - 1) % (ADJUSTED_X_LEN)) + 1
			cur_state = div((transformed - 1) % ((ADJUSTED_X_LEN) * 7), (ADJUSTED_X_LEN)) + 1
			last_codon = div(transformed - 1, ((ADJUSTED_X_LEN) * 7)) + 1
		end

		cur_score = all_scores[first_step]

		if is_x
			cur_score = -Float64(bbt[nuc_start][X_len, aa_pos, cur_state, last_codon])
		else
			cur_score = -Float64(bbt[nuc_start][aa_pos, Y_len, cur_state, last_codon])
		end

		if is_x
			one_path = recurse_back(X_hmm, Y_hmm, X_len, aa_pos, cur_state, cur_score, string(trinucs[last_codon][1]), trinucs[last_codon] * nuc_start, 0.0, [(X_len, aa_pos, cur_state)], L_values, thresh, soft_temp)
			push!(results, one_path)
		else
			one_path = recurse_back(X_hmm, Y_hmm, aa_pos, Y_len, cur_state, cur_score, string(trinucs[last_codon][1]), trinucs[last_codon] * nuc_start, 0.0, [(aa_pos, Y_len, cur_state)], L_values, thresh, soft_temp)
			push!(results, one_path)
		end
	end
	return results
end

function sample_start_point(X_hmm, X_prot, X_len, X_tail_log2d, Y_hmm, Y_prot, Y_len, Y_tail_log2d, num_samples, thresh=1.0, soft_temp=1.0, first_softmax_temp=0.01, gen_bbt=false)
	#Go along all endpoints in all letters and probabilistically choose a start point based on scores.
	vecs_to_sample_x = Dict{AbstractString, Array{Float64, 1}}()
	vecs_to_sample_y = Dict{AbstractString, Array{Float64, 1}}()
	letter_legal_x = Dict{AbstractString, Array{Float64, 1}}()
	letter_legal_y = Dict{AbstractString, Array{Float64, 1}}()

	tail_offset = 0.0

	if !(gen_bbt == false)
		bbt = gen_bbt
	end

	MIN_STOP_OFFSET = 2
	#This is a bit of trickery: our model has a hard time with double stop-codons (because overlapping stop codons are impossible)
	#so we set a minimum offset to keep away from the corner case of X_len, Y_len end points that are impossible to finish.
	#Should be set to 2. That works well. Setting it to 0 or 1 wouldn't work so well (at least in my experience).

	L_values = Dict{Any, Any}()

	for letter_iterator in ["A", "C", "G", "T"]
		top_L = bbt[letter_iterator][1:end, 1:end, 1:end, 1:end]
		for integer in 1:64
			top_L[1:X_len, Y_len:Y_len, 1:7, integer:integer] += reshape(repeat(X_tail_log2d, 1, 7), X_len, 1, 7)
			top_L[X_len:X_len, 1:Y_len, 1:7, integer:integer] += reshape(repeat(Y_tail_log2d, 1, 7), 1, Y_len, 7)
		end
		L_values[letter_iterator] = top_L[1:end, 1:end, 1:end, 1:end]

		final_x_points = top_L[X_len, 1:(Y_len - MIN_STOP_OFFSET), 1:7, 1:64]
		final_y_points = top_L[1:(X_len - MIN_STOP_OFFSET), Y_len, 1:7, 1:64]

		vec_to_sample_x = vec(final_x_points)
		vec_to_sample_y = vec(final_y_points)
		vecs_to_sample_x[letter_iterator] = vec_to_sample_x
		vecs_to_sample_y[letter_iterator] = vec_to_sample_y
	end

	all_scores = Float64[]
	for letter_iterator in ["A", "C", "G", "T"]
		for i in vecs_to_sample_x[letter_iterator]
			push!(all_scores, i)
		end
		for j in vecs_to_sample_y[letter_iterator]
			push!(all_scores, j)
		end
	end

	first_steps = math.softmax(all_scores, num_samples, first_softmax_temp)
	sample_x_len = length(vecs_to_sample_x["A"])
	sample_y_len = length(vecs_to_sample_y["A"])
	per_letter = sample_x_len + sample_y_len

	results = Any[]
	first_step_count = 1
	for first_step in first_steps
		first_step_count += 1

		letter_value = div(first_step - 1, per_letter) + 1
		nuc_start = string("ACGT"[letter_value])
		pos_in = ((first_step - 1) % per_letter) + 1
		is_x = pos_in < sample_x_len
		if is_x
			aa_pos = ((pos_in - 1) % (Y_len-MIN_STOP_OFFSET)) + 1
			cur_state = div((pos_in - 1) % ((Y_len-MIN_STOP_OFFSET) * 7), (Y_len - MIN_STOP_OFFSET)) + 1
			last_codon = div(pos_in - 1, ((Y_len - MIN_STOP_OFFSET) * 7)) + 1
		else
			transformed = mod((first_step - 1), per_letter) + 1 - sample_x_len
			aa_pos = ((transformed - 1) % (X_len - MIN_STOP_OFFSET)) + 1
			cur_state = div((transformed - 1) % ((X_len - MIN_STOP_OFFSET) * 7), (X_len - MIN_STOP_OFFSET)) + 1
			last_codon = div(transformed - 1, ((X_len - MIN_STOP_OFFSET) * 7)) + 1
		end

		cur_score = all_scores[first_step]

		if is_x
			cur_score = -Float64(bbt[nuc_start][X_len, aa_pos, cur_state, last_codon])
		else
			cur_score = -Float64(bbt[nuc_start][aa_pos, Y_len, cur_state, last_codon])
		end

		if is_x
			one_path = recurse_back(X_hmm, Y_hmm, X_len, aa_pos, cur_state, cur_score, string(trinucs[last_codon][1]), trinucs[last_codon] * nuc_start, 0.0, [(X_len, aa_pos, cur_state)], L_values, thresh, soft_temp)
			push!(results, one_path)
		else
			one_path = recurse_back(X_hmm, Y_hmm, aa_pos, Y_len, cur_state, cur_score, string(trinucs[last_codon][1]), trinucs[last_codon] * nuc_start, 0.0, [(aa_pos, Y_len-1, cur_state)], bbt, thresh, soft_temp)
			push!(results, one_path)
		end
	end
	return results
end

#end
