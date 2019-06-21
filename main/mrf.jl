"""
mrf : code associated with dealing with Markov Random Fields, calculating pseudolikelihoods, reading inputs, etc.
		We try to accommodate both Gremlin and CCMPred outputs and we use a parameterization of w1/w2 matrices
		that is easier to code with, so some code to glue the formats together.
"""


module mrf

using JLD, LinearAlgebra

include("types.jl")
include("math.jl")
include("bio_seq.jl")
include("utils.jl")

function debug(x)
	println("DEBUG: $x")
end

function err(x)
	println("ERR: $x")
end

function get_gremlin_consensus(file_name, prot_name) #NOTE: extra dependency (optional_consensus).
	#I maintain some protein translations in a file because they align better with hmmalign or are wild-type versions of proteins.
	#The wild-type thing is relevant b/c you'd presumably be inserting most genes into broader wild-type contexts you'd like long-range interactions to be adapted towards.
	gen_con = true
	if "optional_consensus.fa" in readdir()
		gen_con = false
		reference_seqs = bio_seq.load_fasta("optional_consensus.fa")
		if prot_name in keys(reference_seqs)
			return reference_seqs[prot_name]
		else
			gen_con = true
		end
	end

	if gen_con
		#If I don't have a specific sequence in mind, just generate the consensus sequence (based on w1 matrix only) directly from MRF.
		lookup_nodes = Dict{Int64, Int64}() #based on number of weights, determine number of nodes in MRF.
		for i=1:1500 #assume we don't get proteins longer than this.
			lookup_nodes[((i*(i-1)*21*21)/2) + (i*20)] = i
		end
		w1, w2 = read_weights(file_name)
		nNodes = div(size(w1)[1], 21) #lookup_nodes[length(weights)]
		w1 = reshape(w1, (21, nNodes))'[1:end, 1:20] #this is sorta legacy code to get consensus.
		int_dict = Dict{Int64, String}(1 => "A", 2 => "R", 3 => "N", 4 => "D", 5 => "C", 6 => "Q", 7 => "E", 8 => "G", 9 => "H", 10 => "I", 11 => "L", 12 => "K", 13 => "M", 14 => "F", 15 => "P", 16 => "S", 17 => "T", 18 => "W", 19 => "Y", 20 => "V", 21 => "X")
		prot_con = ""
		for i in 1:nNodes
			prot_con *= int_dict[findmax(w1[i,1:end])[2]]
		end

		return prot_con
	end
end

#Split dense form of weights into w1/w2 variables.
function split_weights(w, nNodes, nStates, nEdges)
	w1 = reshape(w[1:(nStates - 1)*nNodes], nNodes, nStates-1)
	offset = (nStates - 1) * nNodes
	w2 = reshape(w[offset+1:offset+nEdges*nStates^2], nStates, nStates, nEdges)
	return w1, w2
end

function read_weights(jld_file)
	weights = load(replace(jld_file, ".csv"=>".jld")) #["weights"]
	#return weights
	w1 = weights["w1"]
	w2 = weights["w2"]
	return w1, w2
end

function read_csv_weights(csv_file) #Old version just read csv file, jld file is slightly more space efficient for storage.
	#For now reading weights from text file
	in_file = open(csv_file)
	line = readline(in_file)
	weights = Float32[]
	while !eof(in_file)
		push!(weights, parse(Float32, line))
		line = readline(in_file)
	end
	push!(weights, parse(Float32, line))
	close(in_file)
	return weights
end

function transform_weights(w1, w2, nNodes) #transforms strange vector of weights into more interpretable matrix of weights.
	new_w1 = hcat(w1, zeros(nNodes))
	new_w1 = reshape(new_w1', nNodes * 21, 1)
	new_w2 = zeros(Float32, nNodes * 21, nNodes * 21)
	edge_counter = 0
	for edge in collect(utils.ifilter(utils.two_lt, utils.big_vec(nNodes)))
		edge_counter += 1
		ei, ej = edge
		square_mat = w2[1:21, 1:21, edge_counter]
		for mi in 1:21
			for mj in 1:21
				new_w2[((ei - 1) * 21) + mj, ((ej - 1) * 21) + mi] = square_mat[mj, mi]
				new_w2[((ej - 1) * 21) + mi, ((ei - 1) * 21) + mj] = square_mat[mj, mi]
			end
		end
	end
	return new_w1, new_w2
end

function standardize_weights(w1, w2, nNodes)
	#This feels pretty deprecatable.
	new_w1, new_w2 = transform_weights(w1, w2, nNodes)
	return map(Float32, new_w1), map(Float32, new_w2)
end

function init_model(grem_file)
	lookup_nodes = Dict{Int64, Int64}()
	for i=1:1500 #assume we don't get proteins longer than this.
		lookup_nodes[((i*(i-1)*21*21)/2) + (i*20)] = i #Finally don't have to give numNodes. Works for even or odd.
	end
	w1, w2 = read_weights(grem_file)

	nNodes = div(size(w1)[1], 21)

	return types.GremModel(w1, w2, nNodes, 0) #not prot or base_logpot. That comes later.
end

function cpu_assess(prot_mat, w1, w2, nNodes, nProt) #removes cuda dependency at expense of speed.
	my_pv_w2 = prot_mat * w2
	my_ull = my_pv_w2 * prot_mat'
	my_pv_w1 = prot_mat * w1
	numerator = my_pv_w1 + diag(my_ull)
	partitions = reshape((w1' .+ my_pv_w2)', 21, nNodes * nProt)
	logZs = sum(reshape(math.mylogsumexp(partitions'), nNodes, nProt), dims=1)'
	all_scores = numerator - logZs
	return all_scores, my_ull, my_pv_w1, my_pv_w2
end

function psl(prot_seq, w1, w2, verbose = false)
	nNodes = length(prot_seq)
	prot_vec = bio_seq.convert_protein(prot_seq)
	pv_w2 = prot_vec * w2
	pv_w1 = prot_vec * w1
	ull = pv_w1 .+ (pv_w2 * prot_vec')
	ull = ull[1] * -1

	if verbose
		debug("The unnormalized LL is $ull")
	end

	new_idea = w1' + pv_w2
	new_idea = reshape(new_idea, 21, nNodes)
	new_logZ = math.mylogsumexp(new_idea')
	if verbose
		debug("\tThe logZ is $(new_logZ' * ones(nNodes, 1))")
	end
	psi_nll = ull .+ new_logZ' * ones(nNodes, 1)
	return psi_nll[1]
end

function basic_energy_calc(prot_seq, w1, w2)
	nNodes = length(prot_seq)
	prot_vec = bio_seq.convert_protein(prot_seq)
	energy = prot_vec * w1 .+ prot_vec * w2 * prot_vec'
	return energy[1] * -1
end

end
