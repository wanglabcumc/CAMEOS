"""
outparser: CAMEOS output parser
"""

using JLD # Should be placed here too to avoid later crash loading JLD file

module main

include("optimize.jl")
include("types.jl") 

using ArgParse, JLD

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
        "mark_gene"
			help = "1st positional argument, name of gene 'mark'"
			default = "cysJ"
        "deg_gene"
			help = "2nd positional argument, name of gene 'deg'"
			default = "infA"
		"runid"
			help = "3rd positional argument, run id code created by CAMEOS"
            default = "LgolUs7k"
        "--frame"
			help = "optional argument, CAMEOS frame"
			default = "p1"     
	end
	return parse_args(s)
end

function outparse_cameos()
    println("=-= CAMEOS output parser =-= v0.1 - Jan 2021 =-= by LLNL =-=")

    # Parse arguments
   	parsed_args = parse_commandline()
	mark_gene = parsed_args["mark_gene"]
	deg_gene = parsed_args["deg_gene"]
	runid = parsed_args["runid"]
	frame = parsed_args["frame"]

    # Get complete path for input and output files
    subdir = string(mark_gene, "_", deg_gene, "_", frame)
    jld_file = string("saved_pop_", runid, ".jld")
    csv_file = string("summary_", runid, ".csv")
    in_path = string("output/", subdir, "/", jld_file)
    out_path = string("output/", subdir, "/", csv_file)

    # Load JLD file
    print("Loading data from JLD file ", in_path, " ... ")
    variants = load(in_path)["variants"];
    println("OK!")

    # Worte CSV file
    print("Saving data to CSV file ", out_path, " : ")
    parsed = 0
    open(out_path, "w") do io
        println(io, string("full_seq,",
                           mark_gene,"_psls,",
                           mark_gene,"_seq,",
                           deg_gene,"_psls,",
                           deg_gene,"_seq"))
        for var in variants
            println(io, var.full_sequence,",",
                    var.mark_prob,",",
                    var.mark_seq,",", #var.mark_nuc,",",
                    var.deg_prob,",",
                    var.deg_seq) #,",",var.deg_nuc)
            parsed += 1
            if parsed % 200 == 0
                print(".")
            end
        end
    end
    println(" OK!")
    println(parsed, " variants parsed for ", runid)
end
    
@time outparse_cameos()
end
