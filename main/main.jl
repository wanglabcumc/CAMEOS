"""
main : the code you run to double-encode... the "main loop" of optimization.
"""

module main

include("utils.jl")
include("types.jl")
include("lookup.jl")
include("std_setup.jl")
include("math.jl")
include("optimize.jl")
include("mrf.jl")

using ArgParse, Logging, JLD, StatsBase, Distributions, Random, Statistics

function convert_to_saveable(cur_pop)
	new_pop = types.SaveChrome[]

	for indiv in cur_pop
		push!(new_pop, types.SaveChrome(indiv.path, indiv.full_sequence,
										indiv.deg_nuc, indiv.deg_map, indiv.deg_trns, indiv.deg_trne, indiv.deg_d, indiv.deg_skip, indiv.deg_insert,
										indiv.mark_nuc, indiv.mark_map, indiv.mark_trns, indiv.mark_trne, indiv.mark_d, indiv.mark_skip, indiv.mark_insert,
										indiv.deg_prob, indiv.deg_base_E, indiv.deg_seq,
										indiv.mark_prob, indiv.mark_base_E, indiv.mark_seq, indiv.first_weight))
	end

	return new_pop
end

"""
set_up_and_optimize is the code that sets up appropriate information about genes we want to target for double-encoding and
orchestrates the initial HMM double-encoding solution and then iteratively improves on that with greedy pseudolikelihood improvements.
Code also logs information through Logging package.
Genes are often referred to as mark/deg. This is historic from times of thinking of "marker" genes and "designated essential gene".
At this point this distinction is not meaningful but is kept to avoid equally arbitrary x_name, y_name conventions.
Arguments:
	mark_name : gene ID for 'mark' gene. Needed as a key for looking up some values associated with genes in text files.
	deg_name : gene ID for 'deg' gene.
	mark/deg_grem : path to MRF parameter files, usually a JLD/CSV file.
	mark/deg_hmm : path to HMM files, usually .hmm files.
	pop_size : size of population, i.e. number of individual HMM solutions to greedily optimize.
	frame : p1/p2.
	max_iter : maximum number of iterations.
	X_range/Y_range: positions along gene to consider for double-encoding... useful if you want to restrict where double-encoding can occur.
						Ranges are defined in terms of the end of the sequence, so if you have a 70 aa sequence you want to start at position 80,
						set the range to be something like 150:151. Generally some buffer around this value (i.e. 150:160) is useful.
						Don't set, or set to false if you don't care where the genes overlap.
"""
function set_up_and_optimize(log_io, rand_barcode, out_path, mark_name, deg_name, mark_grem, deg_grem, mark_hmm, deg_hmm, pop_size, frame, max_iter, X_range=false, Y_range=false; rand_weights = false, actually_mrf = true)

	@debug("Beginning run.")
	@debug(Libc.strftime(time()))
	@debug(Libc.strftime(time()))

	println("The random barcode on this run is: $rand_barcode")
	@debug("The random barcode on this run is: $rand_barcode")

	do_cull = false #culling reduces number of sequences we optimize over time. Just a trick to save time if you want top-performers only.
	full_even_window = 20
	half_window = div(full_even_window, 2) #window statistics to see if scores still decreasing (sometimes used).
	last_few_sf = Float64[]
	last_few_cs = Float64[]

	#Looks up mean/std. dev of family pseudolikelihoods, pre-computed.
	mu_mark, sig_mark = lookup.mu_sig(deg_name, "psls/")
	mu_deg, sig_deg = lookup.mu_sig(mark_name, "psls/")

	@debug("Stat param vars (mark / deg):\t$(mu_mark)\t$(sig_mark)\t$(mu_deg)\t$(sig_deg)")

	if true #replacing try/catch.
		gen_samples = true #this is true if starting from scratch (general case), false if we have some candidates to optimize ahead of time.
		population = types.Chromosome[]
		local mark_grem_prot, deg_grem_prot
		if !gen_samples
			population = load("trpE_population.jld")["population"][1:100] #[1:100]
			mark_gremodel, deg_gremodel = std_setup.short_set_up(mark_name, deg_name, mark_grem, deg_grem, mark_hmm, deg_hmm, length(population), rand_barcode, frame)
			deg_nNodes, mark_nNodes = deg_gremodel.nNodes, mark_gremodel.nNodes
		else
			#Generate population of sampled hmm starting points.
			if X_range == false || Y_range == false #we require both to be there...
				@debug("Doing standard full set up...")
				mark_gremodel, deg_gremodel, population, mark_grem_prot, deg_grem_prot = std_setup.full_set_up(out_path, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem, pop_size, 1200, 1200, rand_barcode, frame)
			else
				@debug("The x range is $X_range")
				@debug("The y range is $Y_range")
				mark_gremodel, deg_gremodel, population, mark_grem_prot, deg_grem_prot = std_setup.full_sample_set_up(out_path, mark_name, mark_hmm, mark_grem, deg_name, deg_hmm, deg_grem, pop_size, 1200, 1200, rand_barcode, X_range, Y_range, frame)
			end
			deg_nNodes, mark_nNodes = deg_gremodel.nNodes, mark_gremodel.nNodes
		end

		#Look through generated candidates, check their fitness values and their correctness.
		success_chrom, fdp, fmp, fitness_values, indie_deg_maps, indie_mark_maps, deg_ull, deg_pv_w1, deg_pv_w2, deg_prot_mat, mark_ull, mark_pv_w1, mark_pv_w2, mark_prot_mat = optimize.assess_founders(population, mark_gremodel.nNodes, mark_gremodel.w1, mark_gremodel.w2, deg_gremodel.nNodes, deg_gremodel.w1, deg_gremodel.w2)

		@debug("Num success... $(length(success_chrom))")

		if actually_mrf
			founding_fitness_values = fitness_values[1:end]

			mark_base_energy = mrf.basic_energy_calc(mark_grem_prot, mark_gremodel.w1, mark_gremodel.w2)
			deg_base_energy = mrf.basic_energy_calc(deg_grem_prot, deg_gremodel.w1, deg_gremodel.w2)

			#Now we load this info into an ExChrome type object.
			cur_pop = types.ExChrome[]
			for i in 1:length(success_chrom)
				old = success_chrom[i]
				if rand_weights
					new_chrom = types.ExChrome(old.path, old.full_sequence, old.deg_nuc, indie_deg_maps[i], old.deg_trns, old.deg_trne,
															old.deg_d, old.deg_skip, old.deg_insert, old.mark_nuc, indie_mark_maps[i], old.mark_trns, old.mark_trne,
															old.mark_d, old.mark_skip, old.mark_insert, fdp[i], deg_base_energy, deg_prot_mat[i:i, 1:end], utils.aa_vec_to_seq(deg_prot_mat[i:i, 1:end]), deg_ull[i], deg_pv_w1[i], deg_pv_w2[i:i, 1:end], fmp[i], mark_base_energy,
															mark_prot_mat[i:i, 1:end], utils.aa_vec_to_seq(mark_prot_mat[i:i, 1:end]), mark_ull[i], mark_pv_w1[i], mark_pv_w2[i:i, 1:end], rand())
				else
					new_chrom = types.ExChrome(old.path, old.full_sequence, old.deg_nuc, indie_deg_maps[i], old.deg_trns, old.deg_trne,
															old.deg_d, old.deg_skip, old.deg_insert, old.mark_nuc, indie_mark_maps[i], old.mark_trns, old.mark_trne,
															old.mark_d, old.mark_skip, old.mark_insert, fdp[i], deg_base_energy, deg_prot_mat[i:i, 1:end], utils.aa_vec_to_seq(deg_prot_mat[i:i, 1:end]), deg_ull[i], deg_pv_w1[i], deg_pv_w2[i:i, 1:end], fmp[i], mark_base_energy,
															mark_prot_mat[i:i, 1:end], utils.aa_vec_to_seq(mark_prot_mat[i:i, 1:end]), mark_ull[i], mark_pv_w1[i], mark_pv_w2[i:i, 1:end], 0.5)
				end
				push!(cur_pop, new_chrom)
			end

			#There are occasionally differences between HMM/MRF models of the proteins, mostly because HMMs are more flexible with insertions.
			#This generates mappings between positions in either data type.
			mark_hmm_to_grem_map = std_setup.get_explicit_mapping(mark_grem_prot, mark_hmm)
			deg_hmm_to_grem_map = std_setup.get_explicit_mapping(deg_grem_prot, deg_hmm)

			@debug("Hey the explicit mapping for mark_hmm_to_grem is $mark_hmm_to_grem_map")
			@debug("Hey the explicit mapping for deg_hmm_to_grem is $deg_hmm_to_grem_map")
			#So now everything is in an extended chromosome object.

			changed_seq = length(cur_pop)
			no_change_iters = 0
			iter = 0;

			#Let's keep track of optimization history.
			history_matrix = zeros(Float32, round(Int64, pop_size * 1.2), max_iter)

			println("Beginning long-range optimization.")
			@debug("About to start optimization.")
			@debug(Libc.strftime(time()))

			stop_early = false #a condition to evaluate every (few) iteration(s) to see if it's worth continuing a run. When true, we clean up and move to next gene pair.
			not_worth_continuing = false

			#Okay let's set up the energy normal distributions...

			deg_energy_mu, deg_energy_sig = lookup.get_energy_params(deg_name, "energies/")
			mark_energy_mu, mark_energy_sig = lookup.get_energy_params(mark_name, "energies/")

			mark_energy_normal = Normal(mark_energy_mu, mark_energy_sig)
			deg_energy_normal = Normal(deg_energy_mu, deg_energy_sig)

			while (iter == 0 || (iter < max_iter)) #&& !(stop_early) # && minimum((fitness_values)) > 400 && no_change_iters < 100)
				sfv = sort(fitness_values)
				#best_50 = sfv[50]

				if iter % 50 == 0
					println("Step $(iter) of $(max_iter)...")
				end
				@debug("Iteration $iter out of $max_iter")
				@debug("Mean fitness is: $(mean((fitness_values)))")
				@debug("Bottom five are $(sfv[1:5])")
				@debug("... and top five are $(sfv[end-5:end]).")
				@debug("Number of sequences that were changed: $changed_seq")
				flush(log_io)

				#iterated conditional modes, multi-site keep-track = icm_multi_kt. keep track = store values of changes to sequence.
				if do_cull
					if iter < div(max_iter, 4)
						#println("$iter : no cutoff")
						cur_pop, changed_seq = optimize.icm_multi_kt(cur_pop, deg_gremodel, mark_gremodel) #, deg_energy_normal, mark_energy_normal)
					elseif iter >= div(max_iter, 4) && iter < div(max_iter, 2)
						#println("$iter : cutoff is $(sfv[div(length(sfv), 2)])")
						cur_pop, changed_seq = optimize.icm_multi_kt(cur_pop, deg_gremodel, mark_gremodel; score_cutoff = sfv[div(length(sfv), 2)])
					elseif iter >= div(max_iter, 2)
						#println("$iter : cutoff is $(sfv[div(length(sfv), 4)])")
						cur_pop, changed_seq = optimize.icm_multi_kt(cur_pop, deg_gremodel, mark_gremodel; score_cutoff = sfv[div(length(sfv), 4)])
					end
				else
					cur_pop, changed_seq = optimize.icm_multi_kt(cur_pop, deg_gremodel, mark_gremodel)
				end

				fitness_values = optimize.assess_pop(cur_pop)
				summed_fitness = sum(fitness_values) #one of the things we should be watching.

				push!(last_few_sf, summed_fitness)
				push!(last_few_cs, changed_seq)

				indiv_count = 1
				for indiv in fitness_values #history matrix stores values of function being optimized over time.
					history_matrix[indiv_count, iter + 1] = indiv
					indiv_count += 1
				end

				#This was originally just code to test the pseudolikelihoods decreasing correctly. It does.
				#Now it's still sorta somewhat useful to log a few values of sequences in both deg/mark as optimization goes.
				if (mod(iter, 50) == 0 && iter != 0)
					@debug("Second opinion, these should be ~0.")
					for i_count in 1:5
						#for the record this is not a good way to generate random numbers.
						#doesn't go from 1 to end.
						i = i_count #round(Int64, rand() * (length(cur_pop) - 1) + 1)

						bel_deg = cur_pop[i].deg_prob
						bel_mrk = cur_pop[i].mark_prob
						#Sort of assuming that end-1 is because of a "*" at stop codon being added.
						if length(strip(cur_pop[i].deg_seq, '*')) == deg_nNodes && length(strip(cur_pop[i].mark_seq, '*')) == mark_nNodes
							cal_deg = mrf.psl(strip(cur_pop[i].deg_seq, '*'), deg_gremodel.w1, deg_gremodel.w2, false)
							cal_mrk = mrf.psl(strip(cur_pop[i].mark_seq, '*'), mark_gremodel.w1, mark_gremodel.w2, false)
							@debug("DEG: $bel_deg vs. $cal_deg")
							@debug("MRK: $bel_mrk vs. $cal_mrk")
							@debug(abs(bel_deg - cal_deg) < 1e-2)
							@debug(abs(bel_mrk - cal_mrk) < 1e-2)
						else
							@debug("Proteins not correct length!")
							@debug(cur_pop[i].mark_seq)
							@debug(cur_pop[i].deg_seq)
						end
					end
				end
				if changed_seq == 0
					no_change_iters += 1
				else
					no_change_iters = 0
				end
				iter += 1
				@debug("\n")
			end

			#So we're done the iterated conditional modes step. Now we report everything...
			@debug("Done while loop...")
			@debug(Libc.strftime(time()))

			@debug("Let us save the history matrix!")
			save("$out_path/$(mark_name)_$(deg_name)_$frame/opt_his_mat_$(rand_barcode).jld", "hist_mat", history_matrix)

			deg_significance = false
			mark_significance = false #this sees if we get any hits worth keeping at all. If yes we save full population.
			all_cal_deg_scores = Float64[]
			all_cal_mark_scores = Float64[]
			#If no we save just a few (top 10) individuals. These are top-6 overall + top-2 in deg and top-2 in mark.
			out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/all_final_fitness_$(rand_barcode).txt", "w")
			write(out_file, "Ind. #\tMark Score\tDeg Score\tMark Sign\tDeg Sign\n")
			for last_indi in 1:length(cur_pop)
				rep_deg = cur_pop[last_indi].deg_prob
				rep_mrk = cur_pop[last_indi].mark_prob
				if length(strip(cur_pop[last_indi].deg_seq, '*')) == deg_nNodes && length(strip(cur_pop[last_indi].mark_seq, '*')) == mark_nNodes
					cal_deg = abs(mrf.psl(strip(cur_pop[last_indi].deg_seq, '*'), deg_gremodel.w1, deg_gremodel.w2, false))
					cal_mrk = abs(mrf.psl(strip(cur_pop[last_indi].mark_seq, '*'), mark_gremodel.w1, mark_gremodel.w2, false))
					special_deg = "-"
					if cal_deg < mu_deg + sig_deg
						special_deg = "***!"
						deg_significance = true
					elseif cal_deg < mu_deg + 2*sig_deg
						special_deg = "**."
						deg_significance = true
					elseif cal_deg < mu_deg + 3*sig_deg
						special_deg = "*?"
						deg_significance = true
					elseif cal_deg < mu_deg + 4*sig_deg #we don't mind keeping these.
						deg_significance = true
					end
					special_mark = "-"
					if cal_mrk < mu_mark + sig_mark
						special_mark = "***!"
						mark_significance = true
					elseif cal_mrk < mu_mark + 2*sig_mark
						special_mark = "**."
						mark_significance = true
					elseif cal_mrk < mu_mark + 3*sig_mark
						special_mark = "*?"
						mark_significance = true
					elseif cal_mrk < mu_mark + 4*sig_mark #we don't mind keeping these.
						mark_significance = true
					end
					push!(all_cal_mark_scores, cal_mrk)
					push!(all_cal_deg_scores, cal_deg)
					write(out_file, "$(last_indi)\t$(cal_mrk)\t$(cal_deg)\t$special_mark\t$special_deg\n")
				else
					push!(all_cal_mark_scores, Inf)
					push!(all_cal_deg_scores, Inf)
					write(out_file, "$(last_indi)\t-1.0\t-1.0\t/\t/\n")
				end
			end
			close(out_file)

			cal_deg_p = sortperm(all_cal_deg_scores)
			cal_mark_p = sortperm(all_cal_mark_scores)
			cal_both_p = sortperm(fitness_values)

			@debug(Libc.strftime(time()))
			if true #deg_significance && mark_significance
				@debug("Let's save our optimized population...")
				saved_pop = convert_to_saveable(cur_pop)
				save("$out_path/$(mark_name)_$(deg_name)_$frame/saved_pop_$(rand_barcode).jld", "variants", saved_pop)
			else #We want just a subset.
				@debug("We'll save 12 interesting members of the population.")
				selected_pop = ExChrome[]
				for ijk in 1:3
					push!(selected_pop, cur_pop[cal_deg_p[ijk]])
					push!(selected_pop, cur_pop[cal_mark_p[ijk]])
				end
				for ijk in 1:6
					push!(selected_pop, cur_pop[cal_both_p[ijk]])
				end
				save("$out_path/$(mark_name)_$(deg_name)_$frame/top_pop_$(rand_barcode).jld", "variants", selected_pop)
			end

			@debug("And we'll also output our top twelve sequences...")
			out_file = open("$out_path/$(mark_name)_$(deg_name)_$frame/top_twelve_$(rand_barcode).fa", "w")

			for ijk in 1:3
				deg_ind = cal_deg_p[ijk]
				write(out_file, ">TOP_DEG_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
				write(out_file, ">TOP_DEG_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
			end
			for ijk in 1:3
				deg_ind = cal_mark_p[ijk]
				write(out_file, ">TOP_MARK_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
				write(out_file, ">TOP_MARK_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
			end
			for ijk in 1:6
				deg_ind = cal_both_p[ijk]
				write(out_file, ">TOP_GEN_d (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].deg_seq)\n")
				write(out_file, ">TOP_GEN_m (ind $deg_ind ) (fit: $(fitness_values[deg_ind])) (deg: $(all_cal_deg_scores[deg_ind])) (mark: $(all_cal_mark_scores[deg_ind]))\n")
				write(out_file, "$(cur_pop[deg_ind].mark_seq)\n")
			end
			close(out_file)
		end

		@debug("REALLY DONE!")

		@debug(Libc.strftime(time()))

	end
end

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
		"commands"
			help = "positional argument, path to file with tab-delimited commands to run (default: example.txt)"
			default = "example.txt"
		"--num"
			help = "optional, experimentally used for running multiple jobs at once"
			default = "0"
		"--threads"
			help = "optional, experimentally used for running multiple jobs at once"
			default = "1"
	end
	return parse_args(s)
end

function run_file()
	parsed_args = parse_commandline()

	command_file = parsed_args["commands"]
	num = parsed_args["num"]
	threads = parsed_args["threads"]

	run_i = parse(Int64, num)
	NUM_THREADS = parse(Int64, threads)

	if NUM_THREADS > 1
		sleep(rand() * 2) #randomly delay start to de-synchronize starts of runs.
	end

	if isfile(command_file)
		in_file = open(command_file)
		in_read = readlines(in_file)
		close(in_file)

		#This is an additional file that records errors independent of individual log files.
		problem_file = open("problem_runs_$(run_i).txt", "w")

		line_count = 0
		for line in in_read #
			line_count += 1
			if line[1] != '#' && (line_count % NUM_THREADS == run_i)
				try
					run_args = split(line, "\t")
					out_dir, short, long, short_jld, long_jld, short_hmm, long_hmm, pop_size, frame, num_iter = run_args

					rand_barcode = Random.randstring()

					the_out_path = "$out_dir/$(short)_$(long)_$frame/"
					if !(isdir(the_out_path))
						run(`mkdir -p $the_out_path`)
					end

					log_io = open("$out_dir/$(short)_$(long)_$frame/log_$(rand_barcode).txt", "w+")

					logger = SimpleLogger(log_io, Logging.Debug)
					global_logger(logger)

					with_logger(logger) do
						pop_size = parse(Int64, pop_size)
						num_iter = parse(Int64, num_iter)
						set_up_and_optimize(log_io, rand_barcode, out_dir, short, long, short_jld, long_jld, short_hmm, long_hmm, pop_size, frame, num_iter)
					end

					flush(log_io)
					close(log_io)
				catch y
					write(problem_file, "$line $y")
				end

			end
		end
		close(problem_file)
	else
		println("$command_file does not exist. Exiting.")
	end

end

run_file()

end
