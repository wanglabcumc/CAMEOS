"""
std_setup : functions called during standard set-up of optimization. Includes CAMEOS_main() which does most of the work
			computing the tensor with appropriate values. Also functions to deal with input/output presentation of all sequences
			and handling mismatches between HMM/MRF models of sequences.
"""

module std_setup

using Logging, StatsBase, Distributions, Unicode

islower(s) = all(c->islowercase(c) | isspace(c), s);


function debug(x)
	println("DEBUG: $x")
end

function err(x)
	println("ERR: $x")
end

include("types.jl")
#include("math.jl")
include("bio_seq.jl")
include("hmm_logic.jl")
include("dyn_frwd.jl")
include("mrf.jl")
#include("lookup.jl")

function CAMEOS_main(x, x_prot, y, y_prot, frame, extra_out = false)
	start_time = time()

	X_len, X_hmm_state, X_hmm_insert, X_hmm_match = create_HMM(x)
	X_hmm = types.wholeHMM(X_hmm_state, X_hmm_insert, X_hmm_match)

	Y_len, Y_hmm_state, Y_hmm_insert, Y_hmm_match = create_HMM(y)
	Y_hmm = types.wholeHMM(Y_hmm_state, Y_hmm_insert, Y_hmm_match)

	X_len += 2
	Y_len += 2 #for the stop and since we're looking one back...

	mark_name = x_prot
	m_pd, m_nd = load_hmm_seq(mark_name, X_hmm, X_len)
	mark_prot = m_pd[mark_name]
	mark_nuc = m_nd[mark_name]

	deg_name = y_prot
	d_pd, d_nd = load_hmm_seq(deg_name, Y_hmm, Y_len)
	deg_prot = d_pd[deg_name]
	deg_nuc = d_nd[deg_name]

	min_len = min(Y_len, X_len) - 1
	x_init = follow_prot_give_log_2d(Y_hmm, deg_prot, Y_len - min_len)
	y_init = follow_prot_give_log_2d(X_hmm, mark_prot, X_len - min_len)

	all_results = fill_in_tensor(X_hmm, Y_hmm, X_len, Y_len, x_init, y_init, mark_prot, deg_prot, min_len, frame)
	if extra_out
		return all_results, X_hmm, Y_hmm, X_len, Y_len, mark_prot, deg_prot, mark_nuc, deg_nuc
	else
		return all_results
	end
end

function align_consensus(prot_seq, hmm_file)
	hmm_con = get_hmm_consensus(hmm_file)

	alignment = read(pipeline(`printf ">pro_se\n$prot_seq\n"`, `hmmalign --outformat=A2M $hmm_file - `), String)
	alignment = join(split(alignment, '\n')[2:end], "")
	alignment = fix_hmmalign(alignment, hmm_con)

	return alignment
end

function i_j_from_path(path)
	first_elem = path[end]
	last_elem = path[1]
	cur_i = first_elem[1]; cur_j = first_elem[2]
	end_i = last_elem[1]; end_j = last_elem[2]
	return cur_i, cur_j, end_i, end_j
end

#This is logic to format the output of the path through the sequence.
function format_output_map(cur_seq, frame, deg_nuc, mark_nuc, deg_map, mark_map, cur_i, cur_j, end_i, end_j, deg_len, mark_len)
	deg_trns = 1 #Translation start for essential gene.
	deg_trne = 1 #Translation end for essential gene.
	mark_trns = 1 #Same as above.
	mark_trne = 1 #Same...

	if frame == "p1"
		final_mark = ""; final_deg = ""; final_seq = ""
		if mark_map[cur_j] > 1
			final_mark *= lowercase(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1]);
			final_seq *= lowercase(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1])
			final_deg *= "A" * uppercase(cur_seq[2:end]);
			final_seq *= "A" * uppercase(cur_seq[2:end]);
			final_mark *= "A" * uppercase(cur_seq[2:end]) #we always start with A.
			mark_trns = 1
			mark_trne = length(final_seq) - 1
			deg_trns = length(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1]) + 1 #Not sure about the +1...
			deg_trne = length(final_seq) - 1
		elseif deg_map[cur_i] > 1
			final_deg *= lowercase(deg_nuc[1:((deg_map[cur_i - 1]) * 3) - 2] * "At")
			final_seq *= lowercase(deg_nuc[1:((deg_map[cur_i - 1]) * 3) - 2] * "At")
			final_deg *= uppercase(cur_seq);
			final_seq *= uppercase(cur_seq);
			final_mark *= uppercase(cur_seq[2:end])
			deg_trns = 1
			deg_trne = length(final_seq) - 1
			mark_trns = length(deg_nuc[1:((deg_map[cur_i - 1]) * 3)]) + 2 #- 1 + the At stuff...
			mark_trne = length(final_seq) - 1
		else
			#I guess it's possible that both sequences are at 1.
			#I think then we just don't have anything in prefix.
			final_seq *= uppercase(cur_seq)
			final_deg *= uppercase(cur_seq)
			final_mark *= uppercase(cur_seq[2:end])

			deg_trns = 1
			mark_trns = 2
			deg_trne = length(final_seq) - 1
			mark_trne = length(final_seq) - 1
		end
		if (end_j in keys(mark_map)) && mark_map[end_j] < maximum(keys(mark_map)) # div(mark_len, 3)
			final_mark *= lowercase(mark_nuc[(mark_map[end_j] * 3) - 2 : end]);
			final_seq *= lowercase(mark_nuc[(mark_map[end_j] * 3) - 2 : end])
			mark_trne = length(final_seq)
		elseif (end_i in keys(deg_map)) && deg_map[end_i] < maximum(keys(deg_map)) #div(deg_len, 3)
			final_deg *= lowercase(deg_nuc[(deg_map[end_i] * 3) - 1 : end]);
			final_seq *= lowercase(deg_nuc[(deg_map[end_i] * 3) - 1 : end])
			deg_trne = length(final_seq)
		end
	elseif frame == "p2"
		final_mark = ""; final_deg = ""; final_seq = "";
		if mark_map[cur_j] > 1
			final_mark *= lowercase(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1])
			final_seq *= lowercase(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1])
			final_mark *= "A" * uppercase(cur_seq[2:end])
			final_seq *= "A" * uppercase(cur_seq[2:end])
			final_deg *= "A" * uppercase(cur_seq[2:end])
			mark_trns = 1
			mark_trne = length(final_seq) - 1
			deg_trns = length(mark_nuc[1:((mark_map[cur_j - 1]) * 3) - 1]) + 1
			deg_trne = length(final_seq) - 1
		elseif deg_map[cur_i] > 1
			final_deg *= lowercase(deg_nuc[1:((deg_map[cur_i - 1]) * 3)])
			final_seq *= lowercase(deg_nuc[1:((deg_map[cur_i - 1]) * 3)])
			final_deg *= "A" * uppercase(cur_seq[2:end])
			final_seq *= "A" * uppercase(cur_seq[2:end])
			final_mark *= uppercase(cur_seq[2:end])
			deg_trns = 1
			deg_trne = length(final_seq) - 1
			mark_trns = length(deg_nuc[1:((deg_map[cur_i - 1]) * 3)]) + 2
			mark_trne = length(final_seq) - 1
		end
		if (end_j in keys(mark_map)) && mark_map[end_j] < div(mark_len, 3)
			final_mark *= lowercase(mark_nuc[(mark_map[end_j] * 3) - 2:end])
			final_seq *= lowercase(mark_nuc[(mark_map[end_j] * 3) - 2:end])
			mark_trne = length(final_seq)
		elseif (end_i in keys(deg_map)) && deg_map[end_i] < div(deg_len, 3)
			final_deg *= lowercase(deg_nuc[(deg_map[end_i] * 3) - 1:end])
			final_seq *= lowercase(deg_nuc[(deg_map[end_i] * 3) - 1:end])
			deg_trne = length(final_seq)
		end
	end
	return final_deg, final_mark, final_seq, deg_trns, deg_trne, mark_trns, mark_trne
end

function NST_get_explicit_mapping(wt_seq, hmm_file)
	#Returns a mapping from HMM nodes to position in a protein sequence.
	prot_seq = align_consensus(wt_seq, hmm_file)
	mapping = Dict{Int64, Int64}()
	cur_prot = 0
	cur_hmm = 0
	for cha in prot_seq
		if islower(cha)
			cur_prot += 1
		elseif cha == '-'
			cur_hmm += 1
		else
			cur_prot += 1
			cur_hmm += 1
		end
		mapping[cur_hmm] = cur_prot
	end
	mapping[cur_hmm + 1] = cur_prot + 1 #for stop codon
	return mapping
end

function NST_full_general_main(x_prot, x_hmm, y_prot, y_hmm, frame, num_samples, prob_thresh, soft_temp, soft_temp_start, abs_bad, fin_x_nuc="", fin_x_prot="", fin_y_nuc="", fin_y_prot="")
	gen_bbt, X_hmm_obj, Y_hmm_obj, X_len, Y_len, x_prot_seq, y_prot_seq, mark_nuc, deg_nuc = CAMEOS_main(x_hmm, x_prot, y_hmm, y_prot, frame, true)
	println("CAMEOS tensor built")
	cumulative_score = 0.0
	failures = 0
	success_count = 0
	min_len = min(Y_len, X_len) - 1
	deg_map = Dict{Int64, Int64}()
	mark_map = Dict{Int64, Int64}()
	if fin_x_prot != ""
		deg_map = NST_get_explicit_mapping(fin_x_prot, x_hmm)
	else
		deg_map = Dict{Int64, Int64}(i => i for i in 1:X_len)
	end
	if fin_y_prot != ""
		mark_map = NST_get_explicit_mapping(fin_y_prot, y_hmm)
	else
		mark_map = Dict{Int64, Int64}(i => i for i in 1:Y_len)
	end

	all_samples = types.SampleNucs[]
	all_results = Any[]
	#best_result is deterministic run of CAMEOS.
	X_tail = tail_prot_give_log_2d(X_hmm_obj, x_prot_seq, X_len)
	Y_tail = tail_prot_give_log_2d(Y_hmm_obj, y_prot_seq, Y_len)

	best_result = sample_start_point(X_hmm_obj, x_prot_seq, X_len, X_tail, Y_hmm_obj, y_prot_seq, Y_len, Y_tail, 2, 1.0, 1.0, 1e-10, gen_bbt)
	for r in best_result
		push!(all_results, r)
		break #just add one.
	end

	results = sample_start_point(X_hmm_obj, x_prot_seq, X_len, X_tail, Y_hmm_obj, y_prot_seq, Y_len, Y_tail, num_samples, prob_thresh, soft_temp, soft_temp_start, gen_bbt)
	for r in results
		push!(all_results, r)
	end

	for res in all_results
		full_seq, path, cumulative_score, cur_score = res
		cur_i, cur_j, end_i, end_j = i_j_from_path(path)

		final_mark, final_deg, final_seq, mark_trns, mark_trne, deg_trns, deg_trne = format_output_map(full_seq, frame, fin_x_nuc, fin_y_nuc, deg_map, mark_map, cur_i, cur_j, end_i, end_j, X_len * 3, Y_len * 3)

		push!(all_samples, types.SampleNucs(final_seq, deg_trns, deg_trne, mark_trns, mark_trne, path))
	end

	return all_samples
end

function NST_full_general_main_range(x_prot, x_hmm, y_prot, y_hmm, frame, num_samples, prob_thresh, soft_temp, soft_temp_start, abs_bad, x_range, y_range, fin_x_nuc="", fin_x_prot="", fin_y_nuc="", fin_y_prot="")
	#This function only goes through a single range.
	gen_bbt, X_hmm_obj, Y_hmm_obj, X_len, Y_len, x_prot_seq, y_prot_seq, mark_nuc, deg_nuc = CAMEOS_main(x_hmm, x_prot, y_hmm, y_prot, frame, true)
	println("CAMEOS tensor built")
	cumulative_score = 0.0
	failures = 0
	success_count = 0
	min_len = min(Y_len, X_len) - 1
	deg_map = Dict{Int64, Int64}()
	mark_map = Dict{Int64, Int64}()
	if fin_x_prot != ""
		deg_map = NST_get_explicit_mapping(fin_x_prot, x_hmm)
	else
		deg_map = Dict{Int64, Int64}(i => i for i in 1:X_len)
	end
	if fin_y_prot != ""
		mark_map = NST_get_explicit_mapping(fin_y_prot, y_hmm)
	else
		mark_map = Dict{Int64, Int64}(i => i for i in 1:Y_len)
	end

	all_samples = types.SampleNucs[]
	all_results = Any[]

	results = sample_start_point_on_range(X_hmm_obj, x_prot_seq, X_len, Y_hmm_obj, y_prot_seq, Y_len, num_samples, x_range, y_range, prob_thresh, soft_temp, soft_temp_start, gen_bbt)
	for r in results
		push!(all_results, r)
	end

	for res in all_results #[1:5]
		full_seq, path, cumulative_score, cur_score = res
		cur_i, cur_j, end_i, end_j = i_j_from_path(path)
		final_mark, final_deg, final_seq, mark_trns, mark_trne, deg_trns, deg_trne = format_output_map(full_seq, frame, fin_x_nuc, fin_y_nuc, deg_map, mark_map, cur_i, cur_j, end_i, end_j, X_len * 3, Y_len * 3) #ength(fin_x_nuc), length(fin_y_nuc))
		push!(all_samples, types.SampleNucs(final_seq, deg_trns, deg_trne, mark_trns, mark_trne, path))
	end

	return all_samples #, new_all_samples
end

function produce_vectors(trns, trne, skip_item, insert_item)
	tallied_skip = Dict{Int64, Int64}()
	tallied_insert = Dict{Int64, Int64}()

	#Being cautious b/c countmap on an empty list gives an error.
	if length(skip_item) > 0
		tallied_skip = countmap(skip_item) #this looks like: {210 => 4, 228 => 2} kind of thing..
	end
	if length(insert_item) > 0
		tallied_insert = countmap(insert_item) #similarly looks like {195 => 3, 204 => 2}
	end

	ts_keys = sort(collect(keys(tallied_skip)))
	ti_keys = sort(collect(keys(tallied_insert)))

	grem_aa = 0
	bp_counter = 0
	skip_counter = 0
	skip_covered = Int64[]
	insert_covered = Int64[] #trn_info is keep_track[name] in main body

	skip_vec = Int64[]
	insert_vec = Tuple{Int64, Int64}[]
	insert_counter = 0
	bp = trns - 1

	while (bp - insert_counter) <= trne
		bp += 1
		if skip_counter > 0
			push!(skip_vec, bp - insert_counter)
			skip_counter -= 1
    		continue
		end

		if bp_counter % 3 == 0
			grem_aa = div(bp_counter, 3)
		end #So at bp_counter=0, we are at cur_aa = 1. Great.

		if (grem_aa in ts_keys) && !(grem_aa in skip_covered) #then this is going to be skipped
			push!(skip_vec, bp - insert_counter)
			push!(skip_covered, grem_aa) #so we just do this once.
			skip_counter = (3 * tallied_skip[grem_aa]) - 1 #-1 because we're currently adding our base.
			bp_counter += 0 #because we're finishing off the codon.
			continue #don't want to evaluate anything else. Don't want bp_counter to increase
		end

		if (grem_aa in ti_keys) && !(grem_aa in insert_covered) #then we're going to add one or more "..." here
			push!(insert_covered, grem_aa)
			push!(insert_vec, (bp - 1 - insert_counter, tallied_insert[grem_aa]))
			insert_counter += (3 * tallied_insert[grem_aa])
		end

		bp_counter += 1
	end

	return skip_vec, insert_vec
end

function try_run_on_SampleNucs(sample_prots, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem)
	#Same code as run_on_SampleNucs with Exception Handling to avoid cases where problems with protein translation lead to a bug.
	#Bugs are rare and tend to occur with nonsense proteins that confuse alignment via hmmalign and result in very weird sequences.
	#We don't care about these anyway, because they'll be low-scoring. We are fine throwing them away.
	mark_grem_prot = mrf.get_gremlin_consensus(mark_grem, mark_name)
	deg_grem_prot = mrf.get_gremlin_consensus(deg_grem, deg_name)
	sample_count = 1
	mark_sub_prots = Dict{Int64, AbstractString}()
	deg_sub_prots = Dict{Int64, AbstractString}()
	keep_track = Dict{Int64, Any}()
	titles = Int64[]
	deg_prot_trans = ""
	mark_prot_trans = ""
	for sample_nuc in sample_prots
		sample_name = sample_count
		try
			deg_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.deg_trns:sample_nuc.deg_trne])
			mark_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.mark_trns:sample_nuc.mark_trne])
		catch Exception
			println("Exception hit")
			continue
		end

		#If we have no problem translating...
		mark_sub_prots[sample_name] = mark_prot_trans
		deg_sub_prots[sample_name] = deg_prot_trans
		keep_track[sample_name] = sample_nuc
		push!(titles, sample_name)
		sample_count += 1
	end
	deg_skip, deg_insert = redo_alignment(align_grem(deg_grem_prot, deg_hmm), sample_prots, false)
	mark_skip, mark_insert = redo_alignment(align_grem(mark_grem_prot, mark_hmm), sample_prots, true)
	ex_sam = types.ExtendedSampleNucs[]
	for name in titles
		trn_info = keep_track[name]
		try
			deg_skip_vec, deg_insert_vec = produce_vectors(trn_info.deg_trns, trn_info.deg_trne, deg_skip[name], deg_insert[name])
			mark_skip_vec, mark_insert_vec = produce_vectors(trn_info.mark_trns, trn_info.mark_trne, mark_skip[name], mark_insert[name])
			new_sn = types.ExtendedSampleNucs(trn_info.final_seq, trn_info.deg_trns, trn_info.deg_trne, "+", trn_info.mark_trns, trn_info.mark_trne, "+", trn_info.path, deg_skip_vec, deg_insert_vec, mark_skip_vec, mark_insert_vec)
			push!(ex_sam, new_sn)
		catch ne_err
			#I don't care about an exception, just won't use that sample.
		end
	end
	return ex_sam, mark_grem_prot, deg_grem_prot
end

#This is a helper function to deal with insertions and other sequence features that we "censor" meaning that
#model isn't aware that we're making slight modifications to sequence with understanding that it won't affect downstream outputs.
function interpret_censorship(nucs, first_skip, first_insert, first_tran, second_skip, second_insert, second_tran)
	global_skip = union(first_skip, second_skip)

	first_map_local_skip = Int64[]
	for mapped_x in first_skip
		push!(first_map_local_skip, mapped_x)
	end
	first_map_global_skip = Int64[]
	for mapped_x in global_skip
		push!(first_map_global_skip, mapped_x)
	end
	mapped_first_insert = Tuple[]
	for (mapped_x, y) in first_insert
		push!(mapped_first_insert, (mapped_x, y))
	end

	if first_tran[3] == "+"
		first_nuc = types.MRF_nuc(nucs[first_tran[1]:first_tran[2]], first_map_global_skip, first_map_local_skip, mapped_first_insert)
	else
		first_nuc = types.MRF_nuc(bio_seq.rev_comp(nucs[first_tran[2]:first_tran[1]]), first_map_global_skip, first_map_local_skip, mapped_first_insert)
	end

	aa_count = 0
	first_map = Dict{Int64, Int64}()

	for c in first_tran[1]:first_tran[2]
		if !(c in first_map_local_skip)
			aa_count += 1
			first_map[c] = div(aa_count-1, 3) + 1
		end
	end

	second_map_local_skip = Int64[]
	for mapped_x in second_skip
		push!(second_map_local_skip, mapped_x)
	end
	second_map_global_skip = Int64[]
	for mapped_x in global_skip
		push!(second_map_global_skip, mapped_x)
	end
	mapped_second_insert = Tuple[]
	for (mapped_x, y) in second_insert
		push!(mapped_second_insert, (mapped_x, y))
	end

	if second_tran[3] == "+"
		second_nuc = types.MRF_nuc(nucs[second_tran[1]:second_tran[2]], second_map_global_skip, second_map_local_skip, mapped_second_insert)
	elseif second_tran[3] == "-"
		second_nuc = types.MRF_nuc(bio_seq.rev_comp(nucs[second_tran[2]:second_tran[1]]), second_map_global_skip, second_map_local_skip, mapped_second_insert)
	end

	second_map = Dict{Int64, Int64}()
	aa_count = 0
	for c in second_tran[1]:second_tran[2]
		if !(c in second_map_local_skip)
			aa_count += 1
			second_map[c] = div(aa_count-1, 3) + 1
		end
	end

	return first_nuc, first_map, second_nuc, second_map
end

function hmm_score_helper(deg_in_read)
	true_lines = Any[]
	for line in deg_in_read
		if line[1] == '#'
			continue
		else
			push!(true_lines, split(line, " ", keepempty=false))
		end
	end
	deg_hmm_scores = Dict{AbstractString, Float64}()
	for x in true_lines
		prot_id = x[3]
		score = parse(Float64, x[6])
		deg_hmm_scores[prot_id] = score
	end
	return deg_hmm_scores
end

function full_set_up(out_path, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem, pop_size, num_samples, bad_threshold, rand_barcode, frame="p1")
	my_prots = NST_full_general_main(mark_name, mark_hmm, deg_name, deg_hmm, frame, pop_size * 50, 0.90, 1.0, 1.0, 1200, gen_hmm_trace(mark_name, mark_hmm)..., gen_hmm_trace(deg_name, deg_hmm)...)

	println("Evaluating HMM seeds")
	@debug("We've generated samples.")

	sample_count = 1
	mark_sub_prots = Dict{Int64, String}()
	deg_sub_prots = Dict{Int64, String}()
	keep_track = Dict{Int64, Any}()
	titles = Int64[]
	for sample_nuc in my_prots
		sample_name = sample_count
		deg_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.deg_trns:sample_nuc.deg_trne], 1)
		mark_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.mark_trns:sample_nuc.mark_trne], 1)
		mark_sub_prots[sample_name] = mark_prot_trans
		deg_sub_prots[sample_name] = deg_prot_trans
		keep_track[sample_name] = sample_nuc
		push!(titles, sample_name)
		sample_count += 1
	end

	deg_lines = String[]
	mark_lines = String[]
	for sam in titles
		if length(strip(deg_sub_prots[sam])) > 1 && length(strip(mark_sub_prots[sam])) > 1 #hmmscan doesn't like empty lines.
			push!(deg_lines, ">sample_$(sam)\n$(deg_sub_prots[sam])\n")
			push!(mark_lines, ">sample_$(sam)\n$(mark_sub_prots[sam])\n")
		end
	end

	deg_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_out_$(rand_barcode).fa", "w")
	write(deg_out_file, "$(join(deg_lines))")
	close(deg_out_file)

	mark_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_out_$(rand_barcode).fa", "w")
	write(mark_out_file, "$(join(mark_lines))")
	close(mark_out_file)

	deg_hmm_run = read(pipeline(`hmmscan --cpu 1 --tblout=$out_path/$(mark_name)_$(deg_name)_$frame/deg_$(rand_barcode).tbl $deg_hmm $out_path/$(mark_name)_$(deg_name)_$frame/deg_out_$(rand_barcode).fa `), String)
	mark_hmm_run= read(pipeline(`hmmscan --cpu 1 --tblout=$out_path/$(mark_name)_$(deg_name)_$frame/mark_$(rand_barcode).tbl $mark_hmm $out_path/$(mark_name)_$(deg_name)_$frame/mark_out_$(rand_barcode).fa `), String)

	#Okay now read the results and see which samples are top quality.
	deg_in_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_$(rand_barcode).tbl")
	deg_in_read = readlines(deg_in_file)
	close(deg_in_file)
	deg_hmm_scores = hmm_score_helper(deg_in_read)

	mark_in_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_$(rand_barcode).tbl")
	mark_in_read = readlines(mark_in_file)
	close(mark_in_file)
	mark_hmm_scores = hmm_score_helper(mark_in_read)

	#Okay, now we just sum up the contributions of each score in deg/mark for total.
	#We will then choose the top-N scoring samples for further work.

	total_sam = Dict{Int64, Float64}()
	for sam in titles
		sam_key = "sample_$sam"
		if sam_key in keys(mark_hmm_scores) && sam_key in keys(deg_hmm_scores)
			total_score = mark_hmm_scores[sam_key] + deg_hmm_scores[sam_key]
			total_sam[sam] = total_score
		end
	end

	total_sam_num = length(collect(keys(total_sam)))

	safety_margin = round(Int64, pop_size * 1.1) #10% more than pop_size in case of random failures.
	#We sort the dictionary and only keep some of the prots.
	sorted_by_score = sort(collect(total_sam), by= x -> -x[2]) #negative so that top is first.
	elite_candidates = [sorted_by_score[i][1] for i in 1:5]

	#top_quart = sorted_by_score[1:div(total_sam_num, 4)]
	#top_half = sorted_by_score[1:div(total_sam_num, 2)]
	top_third = sorted_by_score[1:div(total_sam_num, 3)]

	top_key_pairs = sample(top_third, safety_margin, replace=false) #clearly safety_margin can never be bigger than top_half/top_quart.
	top_candidates = [top_key_pairs[i][1] for i in 1:safety_margin]

	for elite in elite_candidates
		push!(top_candidates, elite)
	end

	top_candidates = reverse(top_candidates) #to get "elites" in first.

	my_top_prots = types.SampleNucs[]
	for topper in top_candidates
		push!(my_top_prots, my_prots[topper]) #should be an int.
	end

	deg_lines = String[]
	mark_lines = String[]

	#Let's output the top samples just to be able to query its diversity.
	for iiiii in top_candidates
		push!(deg_lines, ">sample_$(iiiii) ($(total_sam[iiiii]))\n$(deg_sub_prots[iiiii])\n")
		push!(mark_lines, ">sample_$(iiiii) ($(total_sam[iiiii]))\n$(mark_sub_prots[iiiii])\n")
	end

	deg_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_select_$(rand_barcode).fa", "w")
	write(deg_out_file, "$(join(deg_lines))")
	close(deg_out_file)

	mark_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_select_$(rand_barcode).fa", "w")
	write(mark_out_file, "$(join(mark_lines))")
	close(mark_out_file)

	@debug("From new_sample_traceback, we get $(length(my_top_prots)) samples with which to work.")
	my_ex_prots, mark_grem_prot, deg_grem_prot = try_run_on_SampleNucs(my_top_prots, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem)
	@debug("We set up with $(length(my_ex_prots)) proteins from auto_align_grem.")

	mark_grem_model = mrf.init_model(mark_grem)
	deg_grem_model = mrf.init_model(deg_grem)

	chromosomes = types.Chromosome[]

	iss = Any[]
	deg_prot_seqs = String[]
	mark_prot_seqs = String[]

	successes = 0
	for ex in my_ex_prots
		if successes >= pop_size #this may happen if our sampling from CAMEOS is out of step with our desired population size.
			break
		end
		#try

		deg_nuc, deg_map, mark_nuc, mark_map = interpret_censorship(ex.final_seq, ex.deg_skip, ex.deg_insert, (ex.deg_trns, ex.deg_trne, ex.deg_d), ex.mark_skip, ex.mark_insert, (ex.mark_trns, ex.mark_trne, ex.mark_d))
		deg_trans = uppercase(bio_seq.translate_constrained_maybe_map(deg_nuc, ex.deg_trns, do_map = false))
		mark_trans = uppercase(bio_seq.translate_constrained_maybe_map(mark_nuc, ex.mark_trns, do_map = false))

		if length(strip(deg_trans, '*')) == deg_grem_model.nNodes && length(strip(mark_trans, '*')) == mark_grem_model.nNodes
			push!(deg_prot_seqs, deg_trans)
			push!(mark_prot_seqs, mark_trans)
			push!(iss, (ex.path, ex.final_seq, deg_nuc, deg_map, ex.deg_trns, ex.deg_trne, ex.deg_d, ex.deg_skip, ex.deg_insert, mark_nuc, mark_map, ex.mark_trns, ex.mark_trne, ex.mark_d, ex.mark_skip, ex.mark_insert))
			#I just want to write all these proteins and assess their likelihoods and stuff as a big batch.
			successes += 1
		end
	end

	@debug("We had $successes sampling successes...")

	deg_prot_mat = zeros(Int64, successes, deg_grem_model.nNodes * 21)
	mark_prot_mat = zeros(Int64, successes, mark_grem_model.nNodes * 21)
	for i in 1:successes
		deg_prot_mat[i, 1:end] = bio_seq.convert_protein(deg_prot_seqs[i][1:deg_grem_model.nNodes])
		mark_prot_mat[i, 1:end] = bio_seq.convert_protein(mark_prot_seqs[i][1:mark_grem_model.nNodes])
	end

	deg_probs, deg_ull, deg_pv_w1, deg_pv_w2 = mrf.cpu_assess(deg_prot_mat, deg_grem_model.w1, deg_grem_model.w2, deg_grem_model.nNodes, successes)
	mark_probs, mark_ull, mark_pv_w1, mark_pv_w2 = mrf.cpu_assess(mark_prot_mat, mark_grem_model.w1, mark_grem_model.w2, mark_grem_model.nNodes, successes)

	for i in 1:successes
		new_chr = types.Chromosome(iss[i]..., deg_probs[i], mark_probs[i])
		push!(chromosomes, new_chr)
	end
	return mark_grem_model, deg_grem_model, chromosomes, mark_grem_prot, deg_grem_prot
end

function full_sample_set_up(out_path, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem, pop_size, num_samples, bad_threshold, rand_barcode, mark_range, deg_range, frame="p1")
	my_sampled_prots = NST_full_general_main_range(mark_name, mark_hmm, deg_name, deg_hmm, frame, pop_size * 50, 0.90, 1.0, 1.0, 1200, mark_range, deg_range, gen_hmm_trace(mark_name, mark_hmm)..., gen_hmm_trace(deg_name, deg_hmm)...)

	sample_count = 1
	mark_sub_prots = Dict{String, String}()
	deg_sub_prots = Dict{String, String}()
	keep_track = Dict{String, Any}()
	titles = String[]
	for sample_nuc in my_sampled_prots
		sample_name = "$sample_count"
		deg_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.deg_trns: sample_nuc.deg_trne], 1)
		mark_prot_trans = bio_seq.translate(sample_nuc.final_seq[sample_nuc.mark_trns: sample_nuc.mark_trne], 1)
		mark_sub_prots[sample_name] = mark_prot_trans
		deg_sub_prots[sample_name] = deg_prot_trans
		keep_track[sample_name] = sample_nuc
		push!(titles, sample_name)
		sample_count += 1
	end

	deg_lines = String[]
	mark_lines = String[]
	for sam in titles
		if length(strip(deg_sub_prots[sam])) > 1 && length(strip(mark_sub_prots[sam])) > 1 #hmmscan doesn't like empty lines.
			push!(deg_lines, ">sample_$(sam)\n$(deg_sub_prots[sam])\n")
			push!(mark_lines, ">sample_$(sam)\n$(mark_sub_prots[sam])\n")
		end
	end

	deg_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_out_$(rand_barcode).fa", "w")
	write(deg_out_file, "$(join(deg_lines))")
	close(deg_out_file)

	mark_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_out_$(rand_barcode).fa", "w")
	write(mark_out_file, "$(join(mark_lines))")
	close(mark_out_file)

	deg_hmm_run = read(pipeline(`hmmscan --cpu 1 --tblout=$out_path/$(mark_name)_$(deg_name)_$frame/deg_$(rand_barcode).tbl $deg_hmm $out_path/$(mark_name)_$(deg_name)_$frame/deg_out_$(rand_barcode).fa `), String)
	mark_hmm_run= read(pipeline(`hmmscan --cpu 1 --tblout=$out_path/$(mark_name)_$(deg_name)_$frame/mark_$(rand_barcode).tbl $mark_hmm $out_path/$(mark_name)_$(deg_name)_$frame/mark_out_$(rand_barcode).fa `), String)

	#Okay now read the results and see which samples are top quality.
	deg_in_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_$(rand_barcode).tbl")
	deg_in_read = readlines(deg_in_file)
	close(deg_in_file)
	deg_hmm_scores = hmm_score_helper(deg_in_read)

	mark_in_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_$(rand_barcode).tbl")
	mark_in_read = readlines(mark_in_file)
	close(mark_in_file)
	mark_hmm_scores = hmm_score_helper(mark_in_read)

	total_sam = Dict{String, Float64}()
	for sam in titles
		sam_key = "sample_$sam"
		if sam_key in keys(mark_hmm_scores) && sam_key in keys(deg_hmm_scores)
			total_score = mark_hmm_scores[sam_key] + deg_hmm_scores[sam_key]
			total_sam[sam] = total_score
		end
	end

	total_sam_num = length(collect(keys(total_sam)))

	safety_margin = round(Int64, pop_size * 1.1) #10% more than pop_size in case of random failures.
	#Okay, well now let's sort the dictionary and only keep some of the prots!
	sorted_by_score = sort(collect(total_sam), by= x -> -x[2]) #negative so that top is first.

	num_from_elite = 5

	num_from_general_sample = safety_margin
	num_from_range_sample = 0

	top_third = sorted_by_score[1:div(total_sam_num, 3)]
	top_key_pairs = sample(top_third, num_from_general_sample, replace=false) #clearly safety_margin can never be bigger than top_half/top_quart.
	top_candidates = [top_key_pairs[i][1] for i in 1:num_from_general_sample]

	elite_candidates = [sorted_by_score[i][1] for i in 1:num_from_elite]

	for elite in elite_candidates
		push!(top_candidates, elite)
	end

	top_candidates = reverse(top_candidates) #to get elites in first.

	my_top_prots = types.SampleNucs[]
	for topper in top_candidates
		push!(my_top_prots, keep_track[topper]) #should be an int.
	end

	deg_lines = String[]
	mark_lines = String[]

	#Let's output the top samples just to be able to query its diversity.
	for iiiii in top_candidates
		push!(deg_lines, ">sample_$(iiiii) ($(total_sam[iiiii]))\n$(deg_sub_prots[iiiii])\n")
		push!(mark_lines, ">sample_$(iiiii) ($(total_sam[iiiii]))\n$(mark_sub_prots[iiiii])\n")
	end

	deg_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/deg_select_$(rand_barcode).fa", "w")
	write(deg_out_file, "$(join(deg_lines))")
	close(deg_out_file)

	mark_out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/mark_select_$(rand_barcode).fa", "w")
	write(mark_out_file, "$(join(mark_lines))")
	close(mark_out_file)

	@debug("From new_sample_traceback, we get $(length(my_top_prots)) samples with which to work.")
	my_ex_prots, mark_grem_prot, deg_grem_prot = try_run_on_SampleNucs(my_top_prots, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem)
	@debug("We set up with $(length(my_ex_prots)) proteins from auto_align_grem.")

	mark_grem_model = mrf.init_model(mark_grem)
	deg_grem_model = mrf.init_model(deg_grem) #MAKE SURE THIS IS FLOAT32.

	chromosomes = types.Chromosome[]

	iss = Any[]
	deg_prot_seqs = String[]
	mark_prot_seqs = String[]

	successes = 0
	for ex in my_ex_prots
		if successes >= pop_size #this may happen if our sampling from CAMEOS is out of step with our desired population size.
			break
		end
		deg_nuc, deg_map, mark_nuc, mark_map = interpret_censorship(ex.final_seq, ex.deg_skip, ex.deg_insert, (ex.deg_trns, ex.deg_trne, ex.deg_d), ex.mark_skip, ex.mark_insert, (ex.mark_trns, ex.mark_trne, ex.mark_d))
		deg_trans = uppercase(bio_seq.translate_constrained_maybe_map(deg_nuc, ex.deg_trns, do_map = false))
		mark_trans = uppercase(bio_seq.translate_constrained_maybe_map(mark_nuc, ex.mark_trns, do_map = false))

		if length(strip(deg_trans, '*')) == deg_grem_model.nNodes && length(strip(mark_trans, '*')) == mark_grem_model.nNodes
			push!(deg_prot_seqs, deg_trans)
			push!(mark_prot_seqs, mark_trans)
			push!(iss, (ex.path, ex.final_seq, deg_nuc, deg_map, ex.deg_trns, ex.deg_trne, ex.deg_d, ex.deg_skip, ex.deg_insert, mark_nuc, mark_map, ex.mark_trns, ex.mark_trne, ex.mark_d, ex.mark_skip, ex.mark_insert))
			#I just want to write all these proteins and assess their likelihoods and stuff as a big batch.
			successes += 1
		end
	end

	@debug("We had $successes sampling successes...")

	deg_prot_mat = zeros(Int64, successes, deg_grem_model.nNodes * 21)
	mark_prot_mat = zeros(Int64, successes, mark_grem_model.nNodes * 21)
	for i in 1:successes
		deg_prot_mat[i, 1:end] = bio_seq.convert_protein(deg_prot_seqs[i][1:deg_grem_model.nNodes])
		mark_prot_mat[i, 1:end] = bio_seq.convert_protein(mark_prot_seqs[i][1:mark_grem_model.nNodes])
	end

	deg_probs, deg_ull, deg_pv_w1, deg_pv_w2 = mrf.cpu_assess(deg_prot_mat, deg_grem_model.w1, deg_grem_model.w2, deg_grem_model.nNodes, successes)
	mark_probs, mark_ull, mark_pv_w1, mark_pv_w2 = mrf.cpu_assess(mark_prot_mat, mark_grem_model.w1, mark_grem_model.w2, mark_grem_model.nNodes, successes)

	for i in 1:successes
		#I can just splat the tuple with "..."
		new_chr = types.Chromosome(iss[i]..., deg_probs[i], mark_probs[i])
		push!(chromosomes, new_chr)
	end
	return mark_grem_model, deg_grem_model, chromosomes, mark_grem_prot, deg_grem_prot
end

function short_set_up(mark_name, deg_name, mark_grem, deg_grem, mark_hmm, deg_hmm, pop_size)
	mark_grem_model = mrf.init_model(mark_grem)
	deg_grem_model = mrf.init_model(deg_grem)

	seqs = bio_seq.load_fasta("gremlin_consensus.fa")
	return mark_grem_model, deg_grem_model
end

function get_explicit_mapping(grem_seq, hmm_file)
	prot_seq = align_consensus(grem_seq, hmm_file)
	#prot_seq = join(split(align, "\n")[2:end])
	mapping = Dict{Int64, Int64}()
	cur_grem = 0
	cur_hmm = 0
	for cha in prot_seq
		if islower(cha)
			cur_grem += 1
		elseif cha == '-'
			cur_hmm += 1
		else
			cur_grem += 1
			cur_hmm += 1
		end
		mapping[cur_hmm] = cur_grem
	end
	mapping[cur_hmm + 1] = cur_grem #for stop codon, but we just map to last character.
	mapping[cur_hmm + 2] = cur_grem
	return mapping
end

end
