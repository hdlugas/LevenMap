
using ArgParse
using Glob
using Printf

function print_help()
    println("""
Help Message for create_count_matrix.jl

Usage:
  julia create_count_matrix.jl --input <path>

Arguments:
  --input_dir   Path to directory containing alignment output from all samples (required)
  --help    Show this help message
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
    for req in ["--input_dir"]
        haskey(args, req) || error("Missing required argument: $req. Use --help for usage.")
    end
    return args
end

function main()
    args = parse_command_line()
    input_dir = args["--input_dir"]
    if !endswith(input_dir, "/")
        input_dir *= "/"
    end
    output_txt = input_dir * "count_matrix.txt"

    alignment_files = Glob.glob("*_alignment.txt", input_dir)

    # Create temp files with summed ncRNA counts per sample
    for sample_file in alignment_files
        uniq_file = sample_file * ".uniq.ncRNAs.tmp"
        if isfile(uniq_file)
            rm(uniq_file)
        end

        ncRNA_counts = Dict{String, Int}()
        open(sample_file, "r") do io
            for line in Iterators.drop(eachline(io), 1)
                seq, ncRNA_raw, count_str = split(line)
                count = parse(Int, count_str)
                
                # Split on "&" and assign count to each ncRNA separately
                for ncRNA in split(ncRNA_raw, "&")
                    ncRNA_counts[ncRNA] = get(ncRNA_counts, ncRNA, 0) + count
                end
            end
        end

        open(uniq_file, "w") do io
            for (ncRNA, count) in ncRNA_counts
                println(io, "$ncRNA $count")
            end
        end
    end

    # Get all unique ncRNAs sorted
    ncRNA_set = Set{String}()
    for sample_file in alignment_files
        open(sample_file, "r") do io
            for line in Iterators.drop(eachline(io), 1)
                _, ncRNA_raw, _ = split(line)
                for ncRNA in split(ncRNA_raw, "&")
                    push!(ncRNA_set, ncRNA)
                end
            end
        end
    end
    sorted_ncRNAs = sort(collect(ncRNA_set))

    # Write header
    open(output_txt, "w") do io
        print(io, "ncRNA")
        for sample_file in alignment_files
            sample_name = basename(sample_file)
            sample_name = replace(sample_name, "_alignment.txt" => "")
            print(io, "\t$sample_name")
        end
        println(io)

        # Write count matrix
        for ncRNA in sorted_ncRNAs
            print(io, ncRNA)
            for sample_file in alignment_files
                uniq_file = sample_file * ".uniq.ncRNAs.tmp"
                count = 0
                if isfile(uniq_file)
                    for line in eachline(uniq_file)
                        split_line = split(line)
                        if split_line[1] == ncRNA
                            count = parse(Int, split_line[2])
                            break
                        end
                    end
                end
                print(io, "\t$count")
            end
            println(io)
        end
    end

    for file in Glob.glob("*.uniq.ncRNAs.tmp", input_dir)
        rm(file)
    end
    println("Count matrix written to $output_txt")
end

main()


