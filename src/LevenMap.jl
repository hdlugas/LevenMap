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
  --output           Path to output TXT file (required); alignment details are also written to <output>_details.txt
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


function exit_with_error(msg::String)
    println(stderr, "ERROR: ", msg)
    exit(1)
end


function check_reference_shards(ref_path::String, read_lengths::Vector{Int}, max_dist::Int)
    rebuild_cmd(min_n, max_n, d) = """
      julia BuildReferenceDatabase.jl \\
        --reference <reference_ncRNAs.fasta> \\
        --output $ref_path \\
        --min_nucleotides $min_n \\
        --max_nucleotides $max_n \\
        --max_dist $d"""

    isdir(ref_path) || exit_with_error("Reference database directory not found: $ref_path\nBuild it first with:\n$(rebuild_cmd(16, 120, max_dist))")
    shard_lengths = Dict{Int, Vector{Int}}()
    for file in readdir(ref_path)
        m = match(r"^ref_n_(\d+)_x_(\d+)\.fasta$", file)
        m === nothing && continue
        push!(get!(shard_lengths, parse(Int, m[2]), Int[]), parse(Int, m[1]))
    end
    isempty(shard_lengths) && exit_with_error("No reference shards (i.e. files named ref_n_<length>_x_<max_dist>.fasta) found in $ref_path\nThis directory does not appear to be a database built by BuildReferenceDatabase.jl. Build it with:\n$(rebuild_cmd(16, 120, max_dist))")

    if !haskey(shard_lengths, max_dist)
        db_max_dist = maximum(keys(shard_lengths))
        lengths = shard_lengths[db_max_dist]
        exit_with_error("LevenMap was run with --max_dist $max_dist, but the reference database $ref_path was built with --max_dist $db_max_dist (i.e. it only supports --max_dist 0 to $db_max_dist). Either rerun LevenMap with --max_dist <= $db_max_dist or rebuild the database with:\n$(rebuild_cmd(minimum(lengths), maximum(lengths), max_dist))")
    end

    missing_lengths = sort(filter(len -> !isfile(ref_path*"ref_n_$(len)_x_$(max_dist).fasta"), read_lengths))
    if !isempty(missing_lengths)
        db_min, db_max = extrema(shard_lengths[max_dist])
        exit_with_error("The input FASTQ contains reads of length $(join(missing_lengths, ", ")) nucleotides, but the reference database $ref_path only covers read lengths $db_min to $db_max nucleotides for --max_dist $max_dist. Either filter reads to length $db_min-$db_max (e.g. cutadapt --minimum-length=$db_min --maximum-length=$db_max) or rebuild the database with:\n$(rebuild_cmd(min(db_min, minimum(read_lengths)), max(db_max, maximum(read_lengths)), max_dist))")
    end
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

max_dist >= 0 || exit_with_error("--max_dist must be a non-negative integer (got $max_dist)")

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

read_list = collect(read_counts)
unique_read_lengths = unique(map(x -> length(x[1]), read_list))
check_reference_shards(ref_path, unique_read_lengths, max_dist)
all_references = Dict{Int, Vector{Tuple{String,String}}}()
for len_val in unique_read_lengths
    all_references[len_val] = load_fasta_references(ref_path*"ref_n_$(len_val)_x_$(max_dist).fasta")
end



##################################### run alignment #####################################
#=
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
	ref_string = join(sort(unique(best_ref_names)), "&")
        push!(results_per_thread[thread_id], (read, ref_string, count))
    end
end

# Allocate one bucket per *actual* Julia thread
results_per_thread = [Vector{Tuple{String,String,Int}}() for _ in 1:nthreads()]
Threads.@threads for i in 1:length(read_list)
    tid = threadid()
    #@assert 1 ≤ tid ≤ length(results_per_thread)
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
        ref_string = join(sort!(unique!(best_ref_names)), "&")
        push!(results_per_thread[tid], (read, ref_string, count))
    end
end


results_per_thread = [Vector{Tuple{String,String,Int}}() for _ in 1:nthreads()]
Threads.@threads for i in 1:length(read_list)
    tid = threadid()
    (read, count) = read_list[i]
    read_revcomp = revcomp(read)

    current_references = get(all_references, length(read), Tuple{String,String}[])
    isempty(current_references) && continue

    min_dist = typemax(Int)
    best_ref_names = String[]
    for (ref_name, ref_seq) in current_references
        d_fwd = evaluate(Levenshtein(), read, ref_seq)
        d_rev = evaluate(Levenshtein(), read_revcomp, ref_seq)
        min_d = min(d_fwd, d_rev)
        if min_d < min_dist
            min_dist = min_d
            empty!(best_ref_names); push!(best_ref_names, ref_name)
        elseif min_d == min_dist
            push!(best_ref_names, ref_name)
        end
    end

    if min_dist <= max_dist
        ref_string = join(sort!(unique!(best_ref_names)), "&")
        push!(results_per_thread[tid], (read, ref_string, count))
    end
end
=#

# --- run alignment (chunked) ---
N = length(read_list)
C = max(1, nthreads())
chunk_size = ceil(Int, N / C)
chunks = [(s:min(s+chunk_size-1, N)) for s in 1:chunk_size:N]

results_per_chunk = [Vector{Tuple{String,String,Int,Int,String,Int}}() for _ in 1:length(chunks)]

@sync for (chunk_idx, rng) in enumerate(chunks)
    @spawn begin
        localbuf = Vector{Tuple{String,String,Int,Int,String,Int}}()
        for i in rng
            (read, count) = read_list[i]
            read_revcomp = revcomp(read)

            current_references = get(all_references, length(read), Tuple{String,String}[])
            isempty(current_references) && continue

            min_dist = typemax(Int)
            best_ref_names = String[]
            best_ref_strands = String[]
            for (ref_name, ref_seq) in current_references
                d_fwd = evaluate(Levenshtein(), read, ref_seq)
                d_rev = evaluate(Levenshtein(), read_revcomp, ref_seq)
                min_d = min(d_fwd, d_rev)
                strand = d_fwd == d_rev ? "+/-" : (d_fwd < d_rev ? "+" : "-")
                if min_d < min_dist
                    min_dist = min_d
                    empty!(best_ref_names); push!(best_ref_names, ref_name)
                    empty!(best_ref_strands); push!(best_ref_strands, strand)
                elseif min_d == min_dist
                    push!(best_ref_names, ref_name)
                    push!(best_ref_strands, strand)
                end
            end
            if min_dist <= max_dist
                # combine strands of references sharing a name so they stay aligned with the sorted unique names
                strand_of = Dict{String,String}()
                for (name, strand) in zip(best_ref_names, best_ref_strands)
                    strand_of[name] = haskey(strand_of, name) && strand_of[name] != strand ? "+/-" : strand
                end
                ref_names = sort!(collect(keys(strand_of)))
                ref_string = join(ref_names, "&")
                strand_string = join([strand_of[name] for name in ref_names], "&")
                push!(localbuf, (read, ref_string, count, min_dist, strand_string, length(ref_names)))
            end
        end
        results_per_chunk[chunk_idx] = localbuf
    end
end




##################################### write output to disk #####################################
open(output_file, "w") do io
    println(io, "Read\tReference\tCount\tEditDistance\tStrand\tNumOptimalAlignments\tMultipleOptimalAlignments")
    for chunk in results_per_chunk
        for (read, ref, count, dist, strand, n_optimal) in chunk
            println(io, "$read\t$ref\t$count\t$dist\t$strand\t$n_optimal\t$(n_optimal > 1)")
        end
    end
end

