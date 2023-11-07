#!/bin/bash

###############
## FUNCTIONS ##
###############

function subsample {
    local file="${1}"
    local sample_size="${2}"
    local subsampled_out="subsampled_${file##*/}"

    seqtk sample -s100 "${in_file}" "${sample_size}" > "${out_dir}/subsampled_reads/${subsampled_out}"
}

function run_skmer {
    local in_file="${1}"
    local skmer_lib="${2}"

    if [ -d "${out_dir}/skmer_library" ]; then
        skmer query -a "${in_file}" "${skmer_lib}/skmer_library" -p "${threads}" 
    else
      	temp_input=$(mktemp -d "${skmer_lib}/temp_in.XXXXXX")
        ln -s $(realpath "${in_file}") "${temp_input}"

        skmer reference "${temp_input}" -l "${skmer_lib}/skmer_library" -p "${threads}"
        rm -r "${temp_input}"
    fi
}

function run_respect {
    echo "hi"
}

function check_coverage {
    local param_file="${1}"

    for x in $(printf $(cat "${param_file}" | cut -f1,2 | tail -n+1 | awk '$2 >= ${high_cov}')); do
        printf "$(printf "${x}" | sed -e 's/ /,/g') "
    done
    for x in $(printf $(cat "${param_file}" | cut -f1,2 | tail -n+1 | awk '$2 <= ${low_cov}')); do
        printf "$(printf "${x}" | sed -e 's/ /,/g') "
    done
}

function get_read_length {
    local dat_file="${1}"
    local read_len=$(grep "read_length" "${dat_file}" | cut -f2)

    if [ ${read_len} == "NA" ]; then
        printf "-1"
    else
        printf "${read_len}"
    fi
}

############
## INPUTS ##
############

input=$1
out_dir="./OUT_subsample_and_estimate"
threads=2
low_cov=3
high_cov=5

usage="bash ${BASH_SOURCE[0]} -h [input] [-o output directory] [-t threads]

Runs nuclear read processing pipeline on a batch of reads split into two mates in reference to a constructed library:
    
    Positional Arguments:
    input       Path to an input directory 

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

#################
## MAIN SCRIPT ##
#################

mkdir "${out_dir}"
mkdir "${out_dir}/subsampled_reads/"


for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 > 10485760' | cut -f 2 | xargs); do
    ## Initial Subampling for 28 Million Reads (Parallelized) ##
    ((i=i%threads)); ((i++==0)) && wait
    subsample "${file}" 28000000 &
done

for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 <= 10485760' | cut -f 2 | xargs); do
    ## Creating symlinks for small files ##
    ln --symbolic $(realpath "${file}") "${out_dir}/subsampled_reads/"
done

for file in $(ls "${out_dir}/subsampled_reads/"); do
    run_skmer "${out_dir}/subsampled_reads/${file}" "${out_dir}"
done

mkdir "${out_dir}/respect_data/"
printf "Input    read_length" > "${out_dir}/respect_data/hist_info.txt"

for file in $(ls "${out_dir}/skmer_library/"); do
    read_len=$(get_read_length "${file}/${file}.dat")
    printf "${file}.hist    ${read_len}" >> "${out_dir}/respect_data/hist_info.txt"
    ln --symbolic "$(realpath ${file}/${file}.hist)" "${out_dir}/respect_data/"
done

# run_respect "${out_dir}/respect_data/"
coverages=$(check_coverage ${out_dir}/estimated_parameters.txt)

for sample in $coverages; do
done

# if [ $(du -k "${input}" | cut -f1) -gt 10485760 ]; then
#     ## Initial Subampling for 28 Million Reads ##
#     subsampled_out="subsampled_${input##*/}"
#     subsample "${input}" 28000000 > "${out_dir}/subsampled_reads/${subsampled_out}"

#     run_skmer "${out_dir}/subsampled_reads/${subsampled_out}" "${out_dir}" 
#     coverage=$(check_coverage "${out_dir}/skmer_library/${subsampled_out%.*}/${subsampled_out%.*}.dat")

#     if [ "${coverage}" -lt "${low_cov}" ] || [ "${coverage}" -gt "${high_cov}" ]; then
#         new_sample_size=$(bc <<< "4*(28000000/${coverage})")
#         subsample "${input}" "${new_sample_size}" > "${out_dir}/subsampled_reads/${subsampled_out}" 
#         ### NOTE: What if new sample is larger than OG file?
#         run_skmer "${out_dir}/subsampled_reads/${subsampled_out}" "${out_dir}"        
#     fi

# else
#     run_skmer "${input}" "${out_dir}"
#     coverage=$(check_coverage "${out_dir}/skmer_library/${input%.*}/${input%.*}.dat")

#     if [ "${coverage}" -gt "${high_cov}" ]; then
#         rm --recursive "${out_dir}/skmer_library/${input%.*}/"
#         new_sample_size=$(bc <<< "4*(28000000/${coverage})")
#         subsample "${input}" "${new_sample_size}" > "${out_dir}/subsampled_reads/${subsampled_out}"
#         run_skmer   "${out_dir}/subsampled_reads/${subsampled_out}" "${out_dir}" 
#     fi
# fi