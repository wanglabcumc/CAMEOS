using JLD, HDF5, GZip

function sub_read(w2::Array{Float64, 2}, co1::Int64, co2::Int64, cur_rows::Array{Float64, 2})
	w2[1 + co1 * 21: (co1 + 1) * 21, 1 + co2 * 21: (co2 + 1) * 21] = cur_rows[1:end, 1:end]
	w2[1 + co2 * 21: (co2 + 1) * 21, 1 + co1 * 21: (co1 + 1) * 21] = cur_rows[1:end, 1:end]'
end

function read_raw(raw_file, out_file)
	in_file = open(raw_file)
	raw_read = readlines(in_file)
	close(in_file)

	ncol = findfirst(n -> n[1] == '#', raw_read) - 1

	w1 = zeros(ncol, 20)
	count = 0
	for line in raw_read[1:ncol]
		count += 1
		w1[count, :] = map(xx -> parse(Float64, xx), split(strip(line), "\t"))
	end

	w2 = zeros(ncol * 21, ncol * 21)
	co1 = 0
	co2 = 0
	row = 0
	cur_rows = zeros(21, 21)
	for line in raw_read[(ncol + 1):end]
		if line[1] == '#'
			if !(co1 == 0 && co2 == 0)
				sub_read(w2, co1, co2, cur_rows)
			end

			row = 0
			pos = split(strip(line))
			cur_rows *= 0.0 #reset
			co1 = parse(Int64, pos[2])
			co2 = parse(Int64, pos[3])
		else
			row += 1
			cur_rows[row, :] = map(mm -> parse(Float64, mm), split(strip(line), "\t"))
		end
	end

	w1_x = hcat(w1, zeros(ncol, 1))'
	w1_alt = reshape(w1_x, (21 * ncol, 1))

	save(out_file, "w1", Float64.(w1_alt), "w2", w2)

end

function main()
	raw_file = ARGS[end - 1]
	jld_file = ARGS[end]
	if isfile(raw_file)
		if endswith(raw_file, ".raw")
			if endswith(jld_file, ".jld")
				read_raw(raw_file, jld_file)
			else
				println("Second file should end in .jld")
			end
		else
			println("First file should end in .raw")
		end
	else
		println("Incorrect input, check if raw file exists.")
	end
end

main()
