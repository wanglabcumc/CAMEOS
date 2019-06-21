"""
lookup : keeps track of code that's main function is to look up external data.
"""

module lookup

using DelimitedFiles, Statistics

function get_energy_params(gene_name, param_dir="energies/") #NOTE: extra dependency
	deg_energy_file = open("$param_dir" * "energy_$(gene_name).txt")
	deg_energy_read = readlines(deg_energy_file)
	close(deg_energy_file)

	#energies = deg_energy[strip(split(der, '\t')[2]) for der in deg_energy_read]
	energies = map(xx -> parse(Float64, xx), deg_energy_read)
	mean_energy = mean(energies)
	std_energy = std(energies)
	return mean_energy, std_energy
end

function mu_sig(gene_name, param_dir="psls/") #NOTE: extra dependency
	for csv_file in readdir(param_dir)
		if csv_file == "psls_$(gene_name).txt"
			#Then it's a match
			raw_data = readdlm(param_dir * "/" * csv_file)
			raw_values = map(Float64, raw_data[1:end, 1])
			mu = mean(raw_values)
			sig = std(raw_values)
			return mu, sig
		end
	end
	return -1, -1 #dummy vars for not successful.
end

end
