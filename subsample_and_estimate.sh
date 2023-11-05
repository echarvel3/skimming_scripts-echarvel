#!/bin/bash


## FUNCTIONS ##
function subsample {
    local in_file="${1}"
    local sample_size="${2}"

    seqtk sample -s100 "${in_file}" "${sample_size}"
}

run_skmer_query () {
    local in_file="${1}"
    local skmer_lib="${2}"

    skmer query -a "${in_file}" "${skmer_lib}" -p "${threads}"
    }

run_skmer_reference () {
    local in_file="${1}"
    local skmer_lib="${2}"

    temp_input=$(mktemp -d "temp_in.XXXXX")
    ln -s $(realpath "${in_file}") "${temp_input}"

    skmer reference "${temp_input}" -o "${skmer_lib}" -p "${threads}"
    rm "${temp_input}"
    }

run_respect () {
    echo "hi"
}

check_coverage () {
    echo "10x"
}
################

input=$1
out_dir="./"
threads=2
low_cov=3
high_cov=5

usage="bash ${BASH_SOURCE[0]} -h [input] [-o output directory] [-t threads]

Runs nuclear read processing pipeline on a batch of reads split into two mates in reference to a constructed library:
    
    Positional Arguments:
    input       Path to a single merged paired-end read sequencing file. 

    Optional inputs:
    -h          Display this help message and exit.
    -o          Path to directory of piepline's output. 
    -t          Number of threads to be used by all software in this micropipeline (seqtk sample, Skmer, RESPECT).

"

while getopts "ho:t:" opts 
do
    case "${opts}" in
        h) echo "${usage}"; exit;;
        o) out_dir="${OPTARG}"; exit;;
        t) threads="${OPTARG}";;
        [?]) echo "invalid input param"; exit 1;;
    esac
done

if [ $(du -k "${input}" | cut -f1) -gt 1000000 ]; then
    ## Initial Subampling for 29 Million Reads ##
    subsampled_out="$out_dir"
    echo "--subampling ${in_file} to ${sample_size}"
    subsample "${input}" 29000000 > "${out_dir}/${subsampled_out}"

    if [ -d "${out_dir}/skmer_library" ]; then
        run_skmer_query "${subsampled_out}" "${out_dir}/skmer_library" 
    else
        run_skmer_reference "${subsampled_out}" "${out_dir}/skmer_library" 
    fi

    coverage=$(check_coverage)
    
    if [ "${coverage}" -lt "${low_cov}" ] || [ "${coverage}" -gt "${high_cov}" ]; then
        new_sample_size=20000
        subsample "${input}" "${new_sample_size}" > "${out_dir}/${subsampled_out}"

        if [ -d "${out_dir}/skmer_library" ]; then
            run_skmer_query "${subsampled_out}" "${out_dir}/skmer_library" 
        else
            run_skmer_reference "${subsampled_out}" "${out_dir}/skmer_library" 
        fi
    fi
# else
#     run_skmer
#     run_respect
#     coverage=check_coverage

#     if [ "${coverage}" -lt "${low_cov}" ] || [ "${coverage}" -gt "${high_cov}" ]; then
#         subsample
#         run_skmer
#         run_respect
#     fi
fi


# if [ 1 -eq "$(echo "$cov_val < $upper_bound" | bc)" ] && [ "0" -eq "$assembly_counter" ]; then
# 	      echo "Coverage in range"
#               echo "Respect running"
#               echo "Input	read_length
# unclassified-kra_${genome}.hist	$(grep read_length ${lib_dir}/unclassified-kra_${genome}/unclassified-kra_${genome}.dat |cut -f2)" >${folder}/respect/${genome}/info.txt

#               respect --debug -N ${iterations} --threads ${threads} -I ${folder}/respect/${genome}/info.txt -i ${lib_dir}/unclassified-kra_${genome}/unclassified-kra_${genome}.hist --tmp ${folder}/respect/${genome}/output/tmp -o ${folder}/respect/${genome}/output || true
#               echo "Respect done"
#               rm ${folder}/respect/${genome}/info.txt
#       elif [ "0" -eq "$assembly_counter" ]; then
# 	      ratio_val="$( echo "scale=4; $set_bound / $cov_val" | bc)"
#               echo "Coverage not in range, sample to be downsampled by ${ratio_val}"
#               seqtk sample ${folder}/kraken/unclassified-kra_${genome}.fq $ratio_val> ${folder}/respect/${genome}/downsampled_${genome}_${ratio_val}.fq
#               echo "Running respect on downsampled sample"
#               respect -i ${folder}/respect/${genome}/downsampled_${genome}_${ratio_val}.fq -N ${iterations} --debug --tmp ${folder}/respect/${genome}/output/tmp -o ${folder}/respect/${genome}/output || true
#               echo "Respect done"
#       else
# 	      if [ "$extension" = "gz" ]; then
# 		      echo "Running respect on genome assembly"
#               	      respect -i ${folder}/assemblies/${genome}.fna -N ${iterations} --debug --tmp ${folder}/respect/${genome}/output/tmp -o ${folder}/respect/${genome}/output || true
#               	      echo "Respect done"
#               else
# 		      temp=${folder}/${genome}.${extension}
# 		      input=$(dirname $lib_dir)
# 		      echo "Running respect on genome assembly"
#                       respect -i ${temp} -N ${iterations} --debug --tmp ${extra}/respect/${genome}/output/tmp -o ${extra}/respect/${genome}/output || true
#                       echo "Respect done"
#               fi
# fi
