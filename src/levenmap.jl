using ArgParse
using FASTX
using DataStructures
using StringDistances
using Base.Threads



##################################### processing functions #####################################
function print_help()
    println("""
Help Message for levenmap.jl

Usage:
  julia levenmap.jl --input_fastq <path> --reference <path> --output <path> [--max_dist <float>] [--threads <int>]

Arguments:
  --input_fastq      Path to input FASTQ file (required)
  --reference        Path to reference database output from build_reference_database.jl (required)
  --output           Path to output TXT file (required)
  --max_dist         Maximum Levenshtein distance allowed in alignment (optional, Int, default = 1)
  --threads          Number of threads (optional, Int, default = 1)
  --help             Show this help message
""")
end


function parse_args()
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
    get!(args, "--max_dist", "1")
    get!(args, "--threads", "1")
    for req in ["--input_fastq", "--reference", "--output"]
        haskey(args, req) || error("Missing required argument: $req. Use --help for usage.")
    end
    return args
end


function load_fasta_references(fasta_file::String)::Vector{Tuple{String, String}}
    references = Vector{Tuple{String, String}}()
    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            header = String(FASTA.identifier(record))
            seq = String(FASTA.sequence(record))
            push!(references, (header, seq))
        end
    end
    return references
end


function revcomp(seq::String)
    complement = Dict('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C', 'N'=>'N')
    return join(complement[base] for base in reverse(seq))
end



##################################### parse command-line arguments #####################################
args = parse_args()
fastq_file = args["--input_fastq"]
ref_path = args["--reference"]
output_file = args["--output"]
max_dist = parse(Int, args["--max_dist"])
num_threads = parse(Int, args["--threads"])

if !endswith(ref_path,"/")
    ref_path *= "/"
end

if nthreads() < num_threads
    println("Warning: Julia was started with $(nthreads()) threads, cannot increase to $num_threads")
end



##################################### process read information #####################################
read_counts = DefaultDict{String, Int}(0)
open(FASTQ.Reader, fastq_file) do reader
    for record in reader
        seq = String(FASTQ.sequence(record))
        read_counts[seq] += 1
    end
end

results_per_thread = [Vector{Tuple{String,String,Int}}() for i in 1:num_threads]
read_list = collect(read_counts)
unique_read_lengths = unique(map(x -> length(x[1]), read_list))
all_references = Dict{Int, Vector{Tuple{String,String}}}()
for len_val in unique_read_lengths
    all_references[len_val] = load_fasta_references(ref_path*"ref_n_$(len_val)_x_$(max_dist).fasta")
end



##################################### run alignment #####################################
Threads.@threads for i in 1:length(read_list)
    thread_id = threadid()
    (read, count) = read_list[i]
    read_revcomp = revcomp(read)

    current_references = all_references[length(read)]
    min_dist = typemax(Int)
    best_ref_names = String[]
    for (ref_name, ref_seq) in current_references
        d_fwd = evaluate(Levenshtein(), read, ref_seq)
        d_rev = evaluate(Levenshtein(), read_revcomp, ref_seq)
	    min_d = min(d_fwd, d_rev)
	    if min_d < min_dist
	        min_dist = min_d
    	    empty!(best_ref_names)
	        push!(best_ref_names, ref_name)
	    elseif min_d == min_dist
	        push!(best_ref_names, ref_name)
        end
    end

    if min_dist <= max_dist
        # Store results in thread-local vector to avoid race
	    ref_string = join(sort(unique(best_ref_names)), "&")
        push!(results_per_thread[thread_id], (read, ref_string, count))
    end
end



##################################### write output to disk #####################################
open(output_file, "w") do io
    println(io, "Read\tReference\tCount")
    for thread_results in results_per_thread
        for (read, ref, count) in thread_results
            println(io, "$read\t$ref\t$count")
        end
    end
end

