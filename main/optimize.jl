"""
optimize : functions associated with optimization of the sequences via MRF (through our ICM-like maximization step).
"""

module optimize

include("utils.jl")
include("types.jl")
include("math.jl")
include("bio_seq.jl")
include("mrf.jl")

using Logging, Distributions

loc_dna_to_num = Dict{Char, Int64}('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)
mut_base_options = bio_seq.get_base_dict()
precomputed_dicodons = bio_seq.precompute_dicodons()


function compute_partitions(perturbed_sums::Array{Float64, 2}, w1::Array{Float32, 2}, nNodes::Int64, num_rows)
	naive_partitions = perturbed_sums .+ w1'
	many_partitions = sum(reshape(math.mylogsumexp(reshape(naive_partitions', 21, nNodes * num_rows)'), nNodes, num_rows), dims=1)'
	partition_rows = reshape(many_partitions, num_rows, 1)
	return partition_rows
end

#consider all di-states, keep track.
#Basically given positions for mutations, we consider all possible products of mutations and return best options.
#We "keep-track" in the sense that we calculate the score while avoiding unnecessary computation as most factors to pseudolikelihood remain unchanged by local variations.
function ca_ds_kt(wt_vec, nNodes, w1, w2, wt_ull, pv_w1, pv_w2, mut_start, mut_end, possibilities)
	necessary_rows = length(possibilities)

	did_change = 0
	region = utils.define_region(mut_start, mut_end)
	cur_aas = Int64[]
	sub_prot = wt_vec[region][1:end]
	cur_aa = findnext(isone, sub_prot, 1)
	while cur_aa != nothing
		push!(cur_aas, cur_aa)
		cur_aa = findnext(isone, sub_prot, cur_aa + 1)
	end

	muts = region[1] .+ (cur_aas .- 1)
	mod_aa = length(cur_aas)

	take_off = sum(w2[muts, 1:end], dims=1)
	new_pv_w2 = pv_w2 .- take_off

	base_spot = region[1] - 1
	additions = zeros(Float64, necessary_rows, nNodes * 21)
	additional_count = 1

	mapping = Dict{Tuple{Int64, Int64}, Int64}()
	for j in 1:21
		for i in 1:21
			if (i, j) in possibilities #Hopefully this is relatively fast to check... could probably be more efficient.
				mapping[(i, j)] = additional_count #we're going to use this mapping to figure out partition based on (i, j) location.
				additions[additional_count, 1:end] = w2[base_spot + i, 1:end] + w2[base_spot + 21 + j, 1:end]
				additional_count += 1
			end
		end
	end

	perturbed_sums = additions .+ new_pv_w2

	needed_partitions = compute_partitions(perturbed_sums, w1, nNodes, necessary_rows)

	true_ull = new_pv_w2 * wt_vec' .- (wt_vec[region]' * new_pv_w2[region])
	#println(true_ull)
	ull_mat = ones(21) * true_ull[1] * ones(21)'
	all_aa_ull = ull_mat .+ 2*new_pv_w2[1:1, region[1:21]]'
	all_aa_ull = all_aa_ull .+ 2*new_pv_w2[1:1, region[22:42]]
	all_aa_ull = all_aa_ull .+ 2*w2[region[1:21], region[22:42]]

	#Then we are only interested in some of these, so...
	all_aa_ull_rows = Float64[]
	for j in 1:21
		for i in 1:21
			if (i, j) in possibilities
				push!(all_aa_ull_rows, all_aa_ull[i, j])
			end
		end
	end

	w1_ull = ones(21) * ((pv_w1 - w1[region]' * wt_vec[region]) * ones(21)')
	w1_ull = w1_ull .+ w1[region[1:21]]
	w1_ull = w1_ull .+ w1[region[22:42]]'

	w1_ull_rows = Float64[]
	for j in 1:21
		for i in 1:21
			if (i, j) in possibilities
				push!(w1_ull_rows, w1_ull[i, j])
			end
		end
	end

	full_ull_rows = all_aa_ull_rows + w1_ull_rows
	psl_rows = needed_partitions - full_ull_rows

	return full_ull_rows, w1_ull_rows, perturbed_sums, psl_rows, mapping
end

#This is not currently used in CAMEOS but the general form of ca_ds_kt to non-contiguous and possibly more mutations is...
#prot_vec is current protein in one-hot form, aas_to_change is positions of aas to change, aa_possibilities is the possible new aas you can get
#w1, w2, pv_w1, pv_w2, nNodes is w1, w2 (MRF), prot_vec * w1, prot_vec * w2, and nodes in MRF.
#using Combinatorics
#=function general_mut(prot_vec, aas_to_change, aa_possibilities, w1, w2, pv_w1, pv_w2, nNodes)
	#let's assume these are given. Let's further require they be in order for uniqueness sake.
	#sorted_aas = sort(aas_to_change)
	num_poss = length(aa_possibilities)
	num_to_change = length(aas_to_change)

	aa_regions = [((21 * (aa - 1)) + 1: aa * 21) for aa in aas_to_change]
	selected_aa = vcat(aa_regions...)
	active_aa = find(y -> y in selected_aa && prot_vec[y] == 1.0, 1:length(prot_vec))

	take_off = sum(w2[active_aa, 1:end], dims=1)
	new_pv_w2 = pv_w2 .- take_off
	additions = zeros(Float64, num_poss, nNodes * 21) #not sure how to define considered_changes.
	additional_count = 1

	for change_vec in aa_possibilities
		aa_indices = [a_r[1] - 1 for a_r in aa_regions] .+ change_vec #warning that I would've expected a_r[1] - 1 but we got rid of it since didn't match.
		additions[additional_count, 1:end] = sum([w2[aa_i, :] for aa_i in aa_indices], dims=1)[1] #w2[aa_indices[1], :] + w2[aa_indices[2], :]
		additional_count += 1
	end

	jiggle_sums = additions .+ new_pv_w2

	needed_partitions = compute_partitions(jiggle_sums, w1, nNodes, num_poss)

	value = pv_w1 - w1[selected_aa]' * prot_vec[selected_aa]

	aa_coordinates = [[aa_regions[i][aa_possibilities[j][i]] for i in 1:num_to_change] for j in 1:num_poss]

	w1_ull = value .+ [sum([w1[aa_regions[i][aa_possibilities[j][i]]] for i in 1:num_to_change]) for j in 1:num_poss]
	#w1_ull = value .+ [sum([w1[aa_i] for aa_i in aa_coordinates[j]]) for j in 1:num_poss] is a different untested way.
	true_ull = new_pv_w2 * prot_vec' - (prot_vec[selected_aa]' * new_pv_w2[selected_aa])

	outside_w2 = true_ull .+ [2 * sum([new_pv_w2[aa_i] for aa_i in aa_coordinates[j]]) for j in 1:num_poss]
	inner_w2 = 2 * [sum([w2[coord[1], coord[2]] for coord in combinations(aa_coordinates[j], 2)]) for j in 1:num_poss]

	w2_full = outside_w2 + inner_w2

	full_ull_rows = w2_full + w1_ull
	psl_rows = needed_partitions - full_ull_rows

	return psl_rows
end=#

#We avoid mutating in regions with insertions, deletions, or at start/end. This function defines acceptable mutation locations.
function acceptable_muts(site, mutated_nucs, individual, deg_nNodes, mark_nNodes)
	all_mark_aa_pos = Int64[]; all_deg_aa_pos = Int64[]; not_set = false;

	for each_nuc in mutated_nucs
		if each_nuc in individual.deg_skip || each_nuc in individual.mark_skip
			not_set = true
			empty!(all_mark_aa_pos)
			empty!(all_deg_aa_pos)
			break
		end

		if each_nuc in keys(individual.deg_map) && individual.deg_map[each_nuc] < deg_nNodes
			deg_aa_pos = individual.deg_map[each_nuc]
			push!(all_deg_aa_pos, deg_aa_pos)
		elseif each_nuc < individual.deg_trns || each_nuc >= (individual.deg_trne - 6) #-6 because we don't want to mess up last codon and each site hits two.
			push!(all_deg_aa_pos, -1) #-1 will be our indication of "outside bounds"
		else
			not_set = true
			empty!(all_mark_aa_pos)
			empty!(all_deg_aa_pos)
			break
		end

		if each_nuc in keys(individual.mark_map) && individual.mark_map[each_nuc] < mark_nNodes
			mark_aa_pos = individual.mark_map[each_nuc]
			push!(all_mark_aa_pos, mark_aa_pos)
		elseif each_nuc < individual.mark_trns || each_nuc >= (individual.mark_trne - 6)
			push!(all_mark_aa_pos, -1)
		else
			not_set = true
			empty!(all_mark_aa_pos)
			empty!(all_deg_aa_pos)
			break
		end

		if -1 in all_mark_aa_pos || -1 in all_deg_aa_pos #remove if this case is handled.
			not_set = true
		end

	end

	#I don't want to modify the start ATGs so I reject mutations on aa_pos of 1.
	if 1 in all_deg_aa_pos || 1 in all_mark_aa_pos
		not_set = true
	end

	#This is the infA array condition (we don't modify some of the amino acids just by design of array)
	#if 1 in all_deg_aa_pos || 2 in all_deg_aa_pos || 3 in all_deg_aa_pos || 4 in all_deg_aa_pos || 69 in all_deg_aa_pos || 70 in all_deg_aa_pos || 71 in all_deg_aa_pos || 72 in all_deg_aa_pos || 73 in all_deg_aa_pos | 74 in all_deg_aa_pos
		#not_set = true
	#end

	if !(not_set) && ((minimum(all_deg_aa_pos) != (maximum(all_deg_aa_pos) - 1)) || (minimum(all_mark_aa_pos) != (maximum(all_mark_aa_pos) - 1)))
		not_set = true
	end

	return not_set, all_mark_aa_pos, all_deg_aa_pos
end

function icm_multi_kt(population, deg_grem, mark_grem, deg_normal = false, mark_normal = false; use_energy = false, score_cutoff = 0) #this expects Normal objects for logpdf fun.
	#This will allow for than one nucleotide to be changed at once by doing a 2-way aa psl calculation.
	mark_seqs = AbstractString[]
	deg_seqs = AbstractString[]
	mark_psls = Float32[]
	deg_psls = Float32[]
	changed_seq = 0
	iss = Any[]
	deg_w1 = deg_grem.w1
	deg_w2 = deg_grem.w2
	deg_nNodes = deg_grem.nNodes
	mark_w1 = mark_grem.w1
	mark_w2 = mark_grem.w2
	mark_nNodes = mark_grem.nNodes

	the_debug_counter = 0
	for individual in population
		the_debug_counter += 1
		if div(length(individual.deg_vec), 21) != deg_nNodes || div(length(individual.mark_vec), 21) != mark_nNodes
			@debug("Trouble with individual $the_debug_counter at icm_mlti_kt")
			continue
		end

		if score_cutoff != 0 && (abs(individual.deg_prob) + abs(individual.mark_prob)) > score_cutoff #we can choose to only optimize individuals below a score cut-off past a certain point of optimization.
			continue #we skip this guy, the sequence stops changing.
		end

		num_muts = 1
		full_len = length(individual.full_sequence)
		unif_len = DiscreteUniform(1, full_len - 3)
		ind_seq = individual.full_sequence[1:end]
		new_seq = ""
		for mut in 1:num_muts
			site = rand(unif_len)
			mut_len = 3 #floor(Int64, rand() * 3) + 2 #number of nucleotides modified (from 2-4).
			mutated_nucs = site:(site+mut_len)
			all_deg_aa_pos = Int64[]
			all_mark_aa_pos = Int64[] #potentially more than one aa getting hit.
			not_set = false

			not_set, all_mark_aa_pos, all_deg_aa_pos = acceptable_muts(site, mutated_nucs, individual, deg_nNodes, mark_nNodes)
			while not_set
				site = rand(unif_len)
				mut_len = 3 #4 guarantees two aas hit... more generally, floor(Int64, rand() * 3) + 2 migh work.
				mutated_nucs = site:(site+mut_len)
				not_set, all_mark_aa_pos, all_deg_aa_pos = acceptable_muts(site, mutated_nucs, individual, deg_nNodes, mark_nNodes)
			end

			unique_all_mark_aa_pos = sort(unique(all_mark_aa_pos))
			unique_all_deg_aa_pos = sort(unique(all_deg_aa_pos))

			deg_codon = Char[]
			mark_codon = Char[]
			deg_start_pos = -1
			mark_start_pos = -1
			original_nuc = '?'

			for nuc in sort(collect(keys(individual.deg_map)))
				if individual.deg_map[nuc] in all_deg_aa_pos
					push!(deg_codon, uppercase(individual.full_sequence[nuc]))
				end
				if nuc == site
					deg_start_pos = length(deg_codon) #should still work, in order.
					original_nuc = loc_dna_to_num[uppercase(individual.full_sequence[nuc])]
				end
			end
			for nuc in sort(collect(keys(individual.mark_map)))
				if individual.mark_map[nuc] in all_mark_aa_pos
					push!(mark_codon, uppercase(individual.full_sequence[nuc]))
				end
				if nuc == site
					mark_start_pos = length(mark_codon)
				end
			end

			resulting_deg_aas = Tuple{Int64, Int64}[]
			resulting_mark_aas = Tuple{Int64, Int64}[]
			mut_dna_options = mut_base_options[mut_len]

			original_nucs = ""
			dna_nuc_counter = 1
			for dna_nucs in mut_dna_options

				tmp_deg_codon = join(deg_codon[1:(deg_start_pos - 1)]) * (dna_nucs) * join(deg_codon[(deg_start_pos + mut_len):end])
				if tmp_deg_codon == join(deg_codon)
					original_nucs = dna_nuc_counter
				end

				push!(resulting_deg_aas, precomputed_dicodons[tmp_deg_codon])
				dna_nuc_counter += 1
			end
			dna_nuc_counter = 1
			for dna_nucs in mut_dna_options

				tmp_mark_codon = join(mark_codon[1:(mark_start_pos - 1)]) * (dna_nucs) * join(mark_codon[(mark_start_pos + mut_len):end])
				if tmp_mark_codon == join(mark_codon)
					original_nucs = dna_nuc_counter
				end

				push!(resulting_mark_aas, precomputed_dicodons[tmp_mark_codon])
				dna_nuc_counter += 1
			end

			deg_possibilities = filter(z -> !(21 in z), unique(resulting_deg_aas))
			mark_possibilities = filter(z -> !(21 in z), unique(resulting_mark_aas))

			#And then it becomes time to actually compute all these options via some _kt() function.
			deg_ull_opt, deg_pv_w1_opt, deg_pv_w2_opt, deg_psl_opt, deg_cads_map = ca_ds_kt(individual.deg_vec, deg_nNodes, deg_w1, deg_w2, individual.deg_ull, individual.deg_pv_w1, individual.deg_pv_w2, unique_all_deg_aa_pos[1], unique_all_deg_aa_pos[2], deg_possibilities)

			mark_ull_opt, mark_pv_w1_opt, mark_pv_w2_opt, mark_psl_opt, mark_cads_map = ca_ds_kt(individual.mark_vec, mark_nNodes, mark_w1, mark_w2, individual.mark_ull, individual.mark_pv_w1, individual.mark_pv_w2, unique_all_mark_aa_pos[1], unique_all_mark_aa_pos[2], mark_possibilities)

			final_options = zeros(Float32, length(mut_dna_options))
			for iii in 1:length(mut_dna_options)
				if !(21 in resulting_deg_aas[iii]) && !(21 in resulting_mark_aas[iii]) #no stop codons allowed!
					if use_energy && mark_normal && deg_normal
						#We multiply by 2.0 so that the "weight" is more intuitively set from 0-1 but values remain reflective of real psls.
						final_options[iii] = 2.0 * individual.first_weight * abs(logpdf(mark_normal, mark_psl_opt[mark_cads_map[resulting_mark_aas[iii]]])) #abs(logpdf(mark_normal, mark_psl_opt[mark_cads_map[resulting_mark_aas[iii]]]))
						final_options[iii] += 2.0 * (1.0 - individual.first_weight) * abs(logpdf(deg_normal, deg_psl_opt[deg_cads_map[resulting_deg_aas[iii]]]))
					else
						final_options[iii] = 2.0 * individual.first_weight * (mark_psl_opt[mark_cads_map[resulting_mark_aas[iii]]] - individual.mark_base_E)
						final_options[iii] += 2.0 * (1.0 - individual.first_weight) * (deg_psl_opt[deg_cads_map[resulting_deg_aas[iii]]] - individual.deg_base_E)
					end
				else
					final_options[iii] = 1_000_000.0
				end
			end

			min_score, min_choice = findmin(final_options)

			best_deg_aa = resulting_deg_aas[min_choice]
			best_mark_aa = resulting_mark_aas[min_choice]

			new_seq = ind_seq[1:(site - 1)] * mut_dna_options[min_choice] * ind_seq[(site + mut_len):end]

			if new_seq != ind_seq
				changed_seq += 1
			end

			new_deg_dna = new_seq[individual.deg_trns:individual.deg_trne]
			new_mark_dna = new_seq[individual.mark_trns:individual.mark_trne]

			individual.deg_nuc = types.MRF_nuc(new_deg_dna, individual.deg_nuc.skip_sam, individual.deg_nuc.skip_trn, individual.deg_nuc.gap_pos)
			individual.mark_nuc = types.MRF_nuc(new_mark_dna, individual.mark_nuc.skip_sam, individual.mark_nuc.skip_trn, individual.mark_nuc.gap_pos)

			individual.deg_seq = strip(uppercase(bio_seq.translate_constrained_maybe_map(individual.deg_nuc, individual.deg_trns, do_map = false)), '*')
			individual.mark_seq = strip(uppercase(bio_seq.translate_constrained_maybe_map(individual.mark_nuc, individual.mark_trns, do_map = false)), '*')

			individual.deg_vec, individual.deg_ull = bio_seq.convert_protein(individual.deg_seq), deg_ull_opt[deg_cads_map[best_deg_aa]]
			individual.deg_pv_w1, individual.deg_pv_w2 = deg_pv_w1_opt[deg_cads_map[best_deg_aa]], vec(deg_pv_w2_opt[deg_cads_map[best_deg_aa], 1:end])'
			individual.deg_prob = deg_psl_opt[deg_cads_map[best_deg_aa]]

			individual.mark_vec, individual.mark_ull = bio_seq.convert_protein(individual.mark_seq), mark_ull_opt[mark_cads_map[best_mark_aa]]
			individual.mark_pv_w1, individual.mark_pv_w2 = mark_pv_w1_opt[mark_cads_map[best_mark_aa]], vec(mark_pv_w2_opt[mark_cads_map[best_mark_aa], 1:end])'
			individual.mark_prob = mark_psl_opt[mark_cads_map[best_mark_aa]]

			individual.full_sequence = new_seq[1:end]

		end
	end
	return population, changed_seq
end

function mask_prot(prot_vec, region)
	new_prot_vec = prot_vec[1:end]
	new_prot_vec[region] = map(Float64, zeros(region))
	return new_prot_vec'
end

function invert_map(nuc_map)
	inverted_map = Dict{Int64, Array{Int64, 1}}()
	for nuc in sort(collect(keys(nuc_map)))
		resulting_aa = nuc_map[nuc]
		if !(resulting_aa in keys(inverted_map))
			inverted_map[resulting_aa] = Int64[]
			push!(inverted_map[resulting_aa], nuc)
		else
			push!(inverted_map[resulting_aa], nuc)
		end
	end
	return inverted_map
end

#Once we end up with generation of samples from the big tensor, we evaluate their pseudolikelihoods.
#This sets a benchmark and also allows us to compute all the quantities we'll need to conpute these values
#faster in the future.
function assess_founders(population, mark_nNodes, mark_w1, mark_w2, deg_nNodes, deg_w1, deg_w2)
	#Need a specific function for first generation. Hereafter we check partitions in icm_mutate().
	#fitness_values = Dict{Chromosome, Float64}()
	fitness_values = Float32[]
	new_degs = String[]
	new_marks = String[]
	successes = 0
	successful_chromosomes = Any[]
	individual_mark_maps = Dict{Int64, Int64}[]
	individual_deg_maps = Dict{Int64, Int64}[]
	for individual in population
		#Translate the chromosomes so that they are protein sequences...
		#We are now going to explicitly map these proteins in this step.
		new_deg_aa, individual_deg_map = bio_seq.translate_constrained_maybe_map(individual.deg_nuc, individual.deg_trns, do_map = true)
		new_deg_aa = uppercase(new_deg_aa)
		new_mark_aa, individual_mark_map = bio_seq.translate_constrained_maybe_map(individual.mark_nuc, individual.mark_trns, do_map = true)
		new_mark_aa = uppercase(new_mark_aa)

		if length(new_deg_aa) >= deg_nNodes && length(new_mark_aa) >= mark_nNodes
			successes += 1
			push!(new_marks, new_mark_aa)
			push!(new_degs, new_deg_aa)
			push!(successful_chromosomes, individual)
			push!(individual_mark_maps, individual_mark_map)
			push!(individual_deg_maps, individual_deg_map)
		else
			@debug("We do hit an issue!") #don't hit any issues.
		end
	end
	#Fill up a prot_mat with these values...
	mark_prot_mat = zeros(Float32, successes, mark_nNodes * 21)
	deg_prot_mat = zeros(Float32, successes, deg_nNodes * 21)
	for i in 1:successes
		mark_prot_mat[i, 1:end] = bio_seq.convert_protein(new_marks[i])[1:mark_nNodes * 21]
		deg_prot_mat[i, 1:end] = bio_seq.convert_protein(new_degs[i])[1:deg_nNodes * 21]
	end

	#Get their partitions and their fitnesses. We mutate based on this info.
	#full_deg_probs, deg_ull, deg_pv_w1, deg_pv_w2 = cuda_assess(deg_prot_mat, deg_cuda_facts)
	#full_mark_probs, mark_ull, mark_pv_w1, mark_pv_w2 = cuda_assess(mark_prot_mat, mark_cuda_facts)

	full_deg_probs, deg_ull, deg_pv_w1, deg_pv_w2 = mrf.cpu_assess(deg_prot_mat, deg_w1, deg_w2, deg_nNodes, successes)
	full_mark_probs, mark_ull, mark_pv_w1, mark_pv_w2 = mrf.cpu_assess(mark_prot_mat, mark_w1, mark_w2, mark_nNodes, successes)

	for i in 1:successes
		if rand() < 0.05
			@debug("Randomly reporting this... deg is $(full_deg_probs[i]) and mark is $(full_mark_probs[i])")
		end
		push!(fitness_values, -full_deg_probs[i] + -full_mark_probs[i])
	end

	#Also we should sum up each of the probs to be a single fitness value.
	return successful_chromosomes, full_deg_probs, full_mark_probs, fitness_values, individual_deg_maps, individual_mark_maps, deg_ull, deg_pv_w1, deg_pv_w2, deg_prot_mat, mark_ull, mark_pv_w1, mark_pv_w2, mark_prot_mat
end

function assess_pop(population)
	fitness_values = Float32[]
	for i in 1:length(population)
		push!(fitness_values, abs(population[i].deg_prob) + abs(population[i].mark_prob))
	end
	return fitness_values
end

end
