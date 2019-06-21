"""
types : stores types / objects used.
	Mostly just bags of parameters for different models and different sequences with various properties we want to maintain.
"""

module types

mutable struct GremModel #MRF model.
	w1::Array{Float32, 2}
	w2::Array{Float32, 2}
	nNodes::Int
	batch_size::Int
end

mutable struct MRF_nuc #from simulated_annealing_dna_clean.
	nucs::String #no reason not to use letters
	skip_sam::Array{Int64, 1}
	skip_trn::Array{Int64, 1}
	#You should guarantee gaps are in sorted order
	gap_pos::Array{Tuple{Int64, Int64}, 1}
end

mutable struct Chromosome #another type
	path::Array{Tuple{Int64, Int64, Int64}}
	full_sequence::AbstractString
	deg_nuc::MRF_nuc
	deg_map::Dict{Int64, Int64}
	deg_trns::Int64
	deg_trne::Int64
	deg_d::AbstractString #this is direction, +/-
	deg_skip::Array{Int64}
	deg_insert::Array{Tuple{Int64, Int64}}

	mark_nuc::MRF_nuc
	mark_map::Dict{Int64, Int64}
	mark_trns::Int64
	mark_trne::Int64
	mark_d::AbstractString
	mark_skip::Array{Int64}
	mark_insert::Array{Tuple{Int64, Int64}}

	deg_prob::Float32
	mark_prob::Float32
end

mutable struct SampleNucs
	final_seq::AbstractString
	deg_trns::Int64
	deg_trne::Int64
	mark_trns::Int64
	mark_trne::Int64
	path #An array of tuples.
end

mutable struct ExtendedSampleNucs
	final_seq::AbstractString
	deg_trns::Int64
	deg_trne::Int64
	deg_d::AbstractString
	mark_trns::Int64
	mark_trne::Int64
	mark_d::AbstractString
	path #An array of tuples.
	deg_skip
	deg_insert
	mark_skip
	mark_insert
end

mutable struct ExChrome #extended chromosome.
	path::Array{Tuple{Int64, Int64, Int64}}
	full_sequence::AbstractString
	deg_nuc::Any
	deg_map::Dict{Int64, Int64}
	deg_trns::Int64
	deg_trne::Int64
	deg_d::AbstractString
	deg_skip::Array{Int64}
	deg_insert::Array{Tuple{Int64, Int64}}

	mark_nuc::Any
	mark_map::Dict{Int64, Int64}
	mark_trns::Int64
	mark_trne::Int64
	mark_d::AbstractString
	mark_skip::Array{Int64}
	mark_insert::Array{Tuple{Int64, Int64}}

	deg_prob::Float32
	deg_base_E::Float32
	deg_vec
	deg_seq::AbstractString
	deg_ull::Float32
	deg_pv_w1::Float32
	deg_pv_w2::Array{Float32, 2}

	mark_prob::Float32
	mark_base_E::Float32
	mark_vec
	mark_seq::AbstractString
	mark_ull::Float32
	mark_pv_w1::Float32
	mark_pv_w2::Array{Float32, 2}

	first_weight::Float32
end

mutable struct SaveChrome #extended chromosome in saveable form (discard matrices which are easily re-calculable)
	path::Array{Tuple{Int64, Int64, Int64}}
	full_sequence::AbstractString
	deg_nuc::Any #MRF_nuc
	deg_map::Dict{Int64, Int64}
	deg_trns::Int64
	deg_trne::Int64
	deg_d::AbstractString
	deg_skip::Array{Int64}
	deg_insert::Array{Tuple{Int64, Int64}}

	mark_nuc::Any #MRF_nuc
	mark_map::Dict{Int64, Int64}
	mark_trns::Int64
	mark_trne::Int64
	mark_d::AbstractString
	mark_skip::Array{Int64}
	mark_insert::Array{Tuple{Int64, Int64}}

	deg_prob::Float32
	deg_base_E::Float32
	#deg_vec
	deg_seq::AbstractString
	#deg_ull::Float32
	#deg_pv_w1::Float32
	#deg_pv_w2::Array{Float32, 2}

	mark_prob::Float32
	mark_base_E::Float32
	#mark_vec
	mark_seq::AbstractString
	#mark_ull::Float32
	#mark_pv_w1::Float32
	#mark_pv_w2::Array{Float32, 2}

	first_weight::Float32
end

mutable struct wholeHMM
	state_probs::Array{Float64, 2}
	insert_probs::Array{Float64, 2}
	match_probs::Array{Float64, 2}
end

end
