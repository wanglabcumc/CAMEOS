"""
utils: a few miscellaneous helper functions
"""

module utils

#Get positions for proteins between i, j.
function define_region(i, j=0)
	if j == 0
		j = i
	end
	return (21 * (i - 1)) + 1 : 21 * j
end

#Lesser than for two inputs (used to sort through edges).
function two_lt(x::Tuple)
	return x[1] < x[2]
end

function big_vec(n)
	outp = Tuple[]
	for i in 1:n
		for j in 1:n
			push!(outp, (i, j))
		end
	end
	return outp
end

function ifilter(f::Function, itr)
	function _it()
		for i = itr
			if f(i)
				produce(i)
			end
		end
	end
	Task(_it)
end

function top_n(el_tuple, val_keys)
  a_max = findmax(el_tuple)
  max_key::Int = val_keys[a_max[2]]
  max_value = a_max[1]
  return (max_key, max_value)
end

function aa_vec_to_seq(vec)
	#Takes in a wide one-hot vector of protein sequence, converts it back to aa sequence.
	aa_seq = Char[]
	prot_len = div(length(vec), 21)
	for i in 1:prot_len
		sub_vec = vec[((i - 1) * 21) + 1 : i * 21]
		sub_aa = findmax(sub_vec)[2]
		push!(aa_seq, "ARNDCQEGHILKMFPSTWYV-"[sub_aa])
	end
	return join(aa_seq)
end

end
