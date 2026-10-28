
using ArgParse
using Printf

##################################### processing functions #####################################
function print_help()
    println("""
Help Message for build_reference_database.jl

Usage:
  julia build_reference_database.jl --reference <path> --output <path> [--min_nucleotides <int>] [--max_nucleotides <int>] [--max_dist <float>]

Arguments:
  --reference          Path to FASTA file of reference ncRNAs (required)
  --output             Path to output directory containing processed reference database (required)
  --min_nucleotides    Minimum length (in units of nucleotides) of reference sequences to include in reference database (optional, Int, default = 12)
  --max_nucleotides    Maximum length (in units of nucleotides) of reference sequences to include in reference database (optional, Int, default = 30)
  --max_dist           Maximum Levenshtein distance allowed in alignment (optional, Int, default = 5)
  --help               Show this help message
""")
end


function parse_command_line()
    args = Dict{String, String}()
    i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--help" || ARGS[i] == "-h"
            print_help()
            exit(0)
        elseif startswith(ARGS[i], "--")
            if i == length(ARGS)
                error("Missing value for argument: ", ARGS[i])
            end
            args[ARGS[i]] = ARGS[i+1]
            i += 2
        else
            error("Unknown argument: ", ARGS[i])
        end
    end
    get!(args, "--min_nucleotides", "12")
    get!(args, "--max_nucleotides", "30")
    get!(args, "--max_dist", "5")
    for req in ["--reference", "--output"]
        haskey(args, req) || error("Missing required argument: $req. Use --help for usage.")
    end
    return args
end


function ensure_directory(path::String)
    isdir(path) || mkpath(path)
end


function filter_fasta(fasta_file::String, output_file::String, min_len::Int, max_len::Int)
    infile = open(fasta_file, "r")
    outfile = open(output_file, "w")

    header = ""
    seq = ""

    for line in eachline(infile)
        line = strip(line)
        if startswith(line, '>')
            if !isempty(seq)
                if length(seq) ≥ min_len && length(seq) ≤ max_len
                    println(outfile, header)
                    println(outfile, seq)
                end
            end
            header = line
            seq = ""
        else
            seq *= line
        end
    end

    # Write last record
    if !isempty(seq) && length(seq) ≥ min_len && length(seq) ≤ max_len
        println(outfile, header)
        println(outfile, seq)
    end

    close(infile)
    close(outfile)
end


function main()
    args = parse_command_line()

    min_nucleotides = parse(Int,args["--min_nucleotides"])
    max_nucleotides = parse(Int,args["--max_nucleotides"])
    max_mismatches = parse(Int,args["--max_dist"])
    fasta_file = args["--reference"]
    output_dir = args["--output"]

    isfile(fasta_file) || error("Reference FASTA file not found: $fasta_file")
    max_mismatches >= 0 || error("--max_dist must be a non-negative integer (got $max_mismatches)")
    min_nucleotides <= max_nucleotides || error("--min_nucleotides ($min_nucleotides) must be <= --max_nucleotides ($max_nucleotides)")

    ensure_directory(output_dir)

    for n_mismatches in 0:max_mismatches
        for n_tmp in min_nucleotides:max_nucleotides
            min_n = n_tmp - n_mismatches
            max_n = n_tmp + n_mismatches
            output_fasta = joinpath(output_dir, @sprintf("ref_n_%d_x_%d.fasta", n_tmp, n_mismatches))
            filter_fasta(fasta_file, output_fasta, min_n, max_n)
        end
    end
end

main()

