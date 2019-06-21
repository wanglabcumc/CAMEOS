"""
dyn_back is code for the traceback. The mxmy, mxiy, mxdy type code simply does the work of taking a HMM state for both
	x/y (MxMy) is "M"atch in X, "M"atch in Y, with "I"nsert, "D"elete other options, calculating the cost of transitioning
	from this state at this node to the next one and thereby either choosing the best scoring transition if we are beneath
	some threshold or otherwise making the choice stochastically through a softmax decision. These functions are what we
	call in the recursive traceback after our tensors have been constructed.
	The numbers 1-7 are used to index the states, which are MxMy, MxIy, MxDy, etc.
"""
module dyn_back

using StatsBase

#This is same as other softmax code but we'll have a local copy.
function back_softmax(vec, num_samples::Int=1, temperature::Float64=1.0)
	vec = vec ./ temperature
	max_value = maximum(vec)

	log_partition = log.(sum(exp.(vec .- max_value))) .+ max_value

	prob_vec = exp.(vec .- log_partition)
	wv = Weights(prob_vec)
	result = sample(1:length(prob_vec), wv, num_samples)
	return result
end

function mxmy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	mxmy_back = bbt[last_letter][short_node - 1, far_node - 1, 1:7, 1:64]

	mxmy_array = zeros(7 * 64)

	#Similar to other functions, we just look at the appropriate next "cell" in tensor and assess scores with the transition before making this "next stop"
	mxmy_array[1:64] = -vec(mxmy_back[1, 1:end]) .+ (short_mm + far_mm)
	mxmy_array[65:128] = -vec(mxmy_back[2, 1:end]) .+ (short_im + far_mm)
	mxmy_array[129:192] = -vec(mxmy_back[3, 1:end]) .+ (short_dm + far_mm)
	mxmy_array[193:256] = -vec(mxmy_back[4, 1:end]) .+ (short_mm + far_im)
	mxmy_array[257:320] = -vec(mxmy_back[5, 1:end]) .+ (short_dm + far_im)
	mxmy_array[321:384] = -vec(mxmy_back[6, 1:end]) .+ (short_mm + far_dm)
	mxmy_array[385:448] = -vec(mxmy_back[7, 1:end]) .+ (short_im + far_dm)

	if rand() < thresh
		next_stop = findmax(-mxmy_array)[2]
	else
		next_stop = back_softmax(-mxmy_array, 1, temperature)[1]
	end
	#some logic about handling this next stop...
	next_state = div(next_stop-1, 64) + 1
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += mxmy_array[next_stop]
	cur_score = mxmy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function mxiy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	mxiy_back = bbt[last_letter][short_node - 1, far_node, 1:5, 1:64]

	mxiy_array = zeros(5 * 64)

	mxiy_array[1:64] = -vec(mxiy_back[1, 1:end]) .+ (short_mm + far_mm)
	mxiy_array[65:128] = -vec(mxiy_back[2, 1:end]) .+ (short_im + far_mm)
	mxiy_array[129:192] = -vec(mxiy_back[3, 1:end]) .+ (short_dm + far_mm)
	mxiy_array[193:256] = -vec(mxiy_back[4, 1:end]) .+ (short_mm + far_im)
	mxiy_array[257:320] = -vec(mxiy_back[5, 1:end]) .+ (short_dm + far_im)

	if rand() < thresh
		next_stop = findmax(-mxiy_array)[2]
	else
		next_stop = back_softmax(-mxiy_array, 1, temperature)[1]
	end

	next_state = div(next_stop-1, 64) + 1
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += mxiy_array[next_stop]
	cur_score = mxiy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function mxdy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	mxdy_back = bbt[last_letter][short_node - 1, far_node - 1, [1, 2, 3, 6, 7], 1:64]

	mxdy_array = zeros(5 * 64)

	mxdy_array[1:64] = -vec(mxdy_back[1, 1:end]) .+ (short_mm + far_mm)
	mxdy_array[65:128] = -vec(mxdy_back[2, 1:end]) .+ (short_im + far_mm)
	mxdy_array[129:192] = -vec(mxdy_back[3, 1:end]) .+ (short_dm + far_mm)
	mxdy_array[193:256] = -vec(mxdy_back[4, 1:end]) .+ (short_mm + far_dm)
	mxdy_array[257:320] = -vec(mxdy_back[5, 1:end]) .+ (short_im + far_dm)

	if rand() < thresh
		next_stop = findmax(-mxdy_array)[2]
	else
		next_stop = back_softmax(-mxdy_array, 1, temperature)[1]
	end

	next_state = [1, 2, 3, 6, 7][div(next_stop-1, 64) + 1][1]
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += mxdy_array[next_stop]
	cur_score = mxdy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function ixmy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	ixmy_back = bbt[last_letter][short_node, far_node - 1, [1, 2, 4, 6, 7], 1:64]

	ixmy_array = zeros(5 * 64)

	ixmy_array[1:64] = -vec(ixmy_back[1, 1:end]) .+ (short_mm + far_mm)
	ixmy_array[65:128] = -vec(ixmy_back[2, 1:end]) .+ (short_im + far_mm)
	ixmy_array[129:192] = -vec(ixmy_back[3, 1:end]) .+ (short_dm + far_mm)
	ixmy_array[193:256] = -vec(ixmy_back[4, 1:end]) .+ (short_mm + far_dm)
	ixmy_array[257:320] = -vec(ixmy_back[5, 1:end]) .+ (short_im + far_dm)

	if rand() < thresh
		next_stop = findmax(-ixmy_array)[2]
	else
		next_stop = back_softmax(-ixmy_array, 1, temperature)[1]
	end

	next_state = [1, 2, 4, 6, 7][div(next_stop-1, 64) + 1][1]
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += ixmy_array[next_stop]
	cur_score = ixmy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function dxmy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	dxmy_back = bbt[last_letter][short_node - 1, far_node - 1, [1, 3, 4, 5, 6], 1:64]
	dxmy_array = zeros(5 * 64)

	dxmy_array[1:64] = -vec(dxmy_back[1, 1:end]) .+ (short_mm + far_mm)
	dxmy_array[65:128] = -vec(dxmy_back[2, 1:end]) .+ (short_dm + far_mm)
	dxmy_array[129:192] = -vec(dxmy_back[3, 1:end]) .+ (short_mm + far_im)
	dxmy_array[193:256] = -vec(dxmy_back[4, 1:end]) .+ (short_dm + far_im)
	dxmy_array[257:320] = -vec(dxmy_back[5, 1:end]) .+ (short_mm + far_dm)

	if rand() < thresh
		next_stop = findmax(-dxmy_array)[2]
	else
		next_stop = back_softmax(-dxmy_array, 1, temperature)[1]
	end

	next_state = [1, 3, 4, 5, 6][div(next_stop - 1, 64) + 1][1]
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += dxmy_array[next_stop]
	cur_score = dxmy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function ixdy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	ixdy_back = bbt[last_letter][short_node, far_node - 1, [1, 2, 6, 7], 1:64]

	ixdy_array = zeros(4 * 64)

	ixdy_array[1:64] = -vec(ixdy_back[1, 1:end]) .+ (short_mm + far_mm)
	ixdy_array[65:128] = -vec(ixdy_back[2, 1:end]) .+ (short_im + far_mm)
	ixdy_array[129:192] = -vec(ixdy_back[3, 1:end]) .+ (short_mm + far_dm)
	ixdy_array[193:256] = -vec(ixdy_back[4, 1:end]) .+ (short_im + far_dm)

	if rand() < thresh
		next_stop = findmax(-ixdy_array)[2]
	else
		next_stop = back_softmax(-ixdy_array, 1, temperature)[1]
	end

	next_state = [1, 2, 6, 7][div(next_stop - 1, 64) + 1][1]
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]
	cumulative_score += ixdy_array[next_stop]
	cur_score = ixdy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

function dxiy(trinucs, last_letter::AbstractString, short_node::Int64, far_node::Int64, cumulative_score::Float64, cur_score::Float64, passed_values::Array{Float64, 1}, bbt, thresh, temperature)
	short_mm, short_md, short_im, short_dm, short_dd, short_mi, short_ii, far_mm, far_md, far_im, far_dm, far_dd, far_mi, far_ii = passed_values;
	dxiy_back = bbt[last_letter][short_node - 1, far_node, [1, 3, 4, 5], 1:64]

	dxiy_array = zeros(4 * 64)

	dxiy_array[1:64] = -vec(dxiy_back[1, 1:end]) .+ (short_mm + far_mm)
	dxiy_array[65:128] = -vec(dxiy_back[2, 1:end]) .+ (short_dm + far_mm)
	dxiy_array[129:192] = -vec(dxiy_back[3, 1:end]) .+ (short_mm + far_im)
	dxiy_array[193:256] = -vec(dxiy_back[4, 1:end]) .+ (short_dm + far_im)

	if rand() < thresh
		next_stop = findmax(-dxiy_array)[2]
	else
		next_stop = back_softmax(-dxiy_array, 1, temperature)[1]
	end

	next_state = [1, 3, 4, 5][div(next_stop - 1, 64) + 1][1]
	next_letter = div(mod(next_stop-1, 64), 16) + 1
	next_codon = trinucs[mod(next_stop-1, 64) + 1]

	cumulative_score += dxiy_array[next_stop]
	cur_score = dxiy_array[next_stop]
	return next_state, next_letter, next_codon, cur_score, cumulative_score
end

end
