"""
math : keeps track of a few mathy functions.
"""

module math

using StatsBase

function mylogsumexp(b)
	B = maximum(b, dims=2)
	return log.(sum(exp.(broadcast(-, b, B)), dims=2)) + B

end

function abs_tol(older, newer, thresh=0.0002)
	if abs(older - newer) < thresh
		return true
	else
		return false
	end
end

function rel_tol(older, newer, thresh=0.0004)
	fun = abs(1 - (older / newer))
	if fun < thresh
		return true
	else
		return false
	end
end

function softmax(vec, num_samples::Int=1, temperature::Float64=1.0)
	vec = vec ./ temperature
	max_value = maximum(vec)
	log_partition = log.(sum(exp.(vec .- max_value))) .+ max_value

	prob_vec = exp.(vec .- log_partition)

	wv = Weights(prob_vec)
	result = sample(1:length(prob_vec), wv, num_samples)
	return result
end

end
