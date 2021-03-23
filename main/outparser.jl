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
            arg_type = String
			default = "cysJ"
        "deg_gene"
		    help = "2nd positional argument, name of gene 'deg'"
            arg_type = String
			default = "infA"
		"runid"
			help = "3rd positional argument, run id code created by CAMEOS"
            arg_type = String
            default = "LgolUs7k"
        "--frame"
			help = "optional argument, CAMEOS frame"
            arg_type = String
		    default = "p1"     
        "--fasta"
            help = "generate a FASTA file too"
            action = :store_true
        "--just-fullseq"
            help = "the FASTA will only contain the full sequence"
            action = :store_true
	end
	return parse_args(s)
end

function outparse_cameos()
    println("=-= CAMEOS output parser =-= v0.4 - Mar 2021 =-= by LLNL =-=")

    # Parse arguments
   	parsed_args = parse_commandline()
	mark_gene = parsed_args["mark_gene"]
	deg_gene = parsed_args["deg_gene"]
	runid = parsed_args["runid"]
	frame = parsed_args["frame"]
    fasta = parsed_args["fasta"]
    justfullseq = parsed_args["just-fullseq"]

    # Get complete path for input and output files
    subdir = string(mark_gene, "_", deg_gene, "_", frame)
    jld_file = string("saved_pop_", runid, ".jld")
    csv_file = string("summary_", runid, ".csv")
    fa_file = string("variants_", runid, ".fasta")
    in_path = string("output/", subdir, "/", jld_file)
    out_path = string("output/", subdir, "/", csv_file)
    fa_path = string("output/", subdir, "/", fa_file)

    # Load JLD file
    print("Loading data from JLD file ", in_path, " ... ")
    variants = load(in_path)["variants"];
    println("OK!")

    # Save CSV file
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

    # Save fasta file
    if fasta
        print("Saving data to FASTA file ", fa_path, " : ")
        variant = 1
        open(fa_path, "w") do io

            function print2fasta(var, field::String, nucs::Bool=false)
                # Prints header and sequence to fasta file
                println(io, ">CAMEOS run:", runid,
                        " mark_gene:", mark_gene,
                        " deg_gene:", deg_gene,
                        " variant:", variant, "/", parsed,
                        " ", string(field))
                if nucs
                    println(io, getfield(getfield(var, Symbol(field)), :nucs))
                else
                    println(io, getfield(var, Symbol(field)))
                end
                return nothing
            end

            for var in variants
                if justfullseq
                    print2fasta(var, "full_sequence")  # equal to mark_nuc.nucs
                else
                    print2fasta(var, "mark_seq")
                    print2fasta(var, "mark_nuc", true)
                    print2fasta(var, "deg_seq")
                    print2fasta(var, "deg_nuc", true)
                end
                variant += 1
                if variant % 200 == 0
                    print(".")
                end
            end
        end
        println(" OK!")
        println(variant-1, " variants in FASTA file for ", runid)    
    end
end
    
@time outparse_cameos()
end
