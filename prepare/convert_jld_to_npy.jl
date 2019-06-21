using NPZ, JLD

function main()
	in_file = ARGS[end]
	if isfile(in_file)
		if endswith(in_file, ".jld")
			model = load(in_file)

			w1 = model["w1"]
			w2 = model["w2"]

			npzwrite(replace(in_file, ".jld" => "_w1") * ".npy", w1)
			npzwrite(replace(in_file, ".jld" => "_w2") * ".npy", w2)
		else
			println("$in_file needs to be a JLD file, have suffix .jld")
		end
	else
		println("$in_file is not a file")
	end
end

main()
