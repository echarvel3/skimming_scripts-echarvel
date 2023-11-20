#!/bin/bash

###############
## FUNCTIONS ##
###############

function pwait {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
        sleep 1
    done
}

function subsample {
    local in_file="${1}"
    local sample_size="${2}"
    local subsampled_out="subsampled_${in_file##*/}"

    seqtk sample -s100 "${in_file}" "${sample_size}" > "${out_dir}/subsampled_reads/${subsampled_out}"
}

function calc_new_subsample {
    local sample=${1}
    local genome=${sample%.hist,*}
    local cov=${sample#*,}
    local optimal_samp=$(bc <<< "${mean_cov}*${initial_sampling}/${cov}")

    echo $genome $cov $optimal_samp $input $out_dir

    if [ ${sample%%_*} == "subsampled" ]; then
        subsample "${input}/${genome#"subsampled_"}.fastq" "${optimal_samp}" "${out_dir}"
    fi
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

function check_coverage {
    local param_file="${1}"

    for x in $(cat "${param_file}" | cut -f1,4 | tail -n+2 | awk -v high_cov="${high_cov}" '$2 > high_cov' | sed -e 's/\t/,/g'); do
        echo "$(echo "${x}")"
    done
    for x in $(cat "${param_file}" | cut -f1,4 | tail -n+2 | awk -v high_cov="${high_cov}" '$2 < low_cov' | sed -e 's/\t/,/g'); do
        echo "$(echo "${x}")"
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

SCRIPT_DIR="/home/echarvel/skimming_scripts/"

input=$1
out_dir="./OUT_subsample_and_estimate"
threads=10
low_cov=3
high_cov=5
mean_cov=4
initial_sampling=30000000
file_ending1="1.fq.gz"

usage="bash ${BASH_SOURCE[0]} -h [input] [-o output directory] [-t threads]

Runs nuclear read processing pipeline on a batch of reads split into two mates in reference to a constructed library:
    
    Positional Arguments:
    input       Path to an input directory 

    Optional inputs:
    -h          Display this help message and exit.
    -o          Path to directory of piepline's output. 
    -t          Threads to be used by all software in this micropipeline (seqtk sample, Skmer, RESPECT).
    -m          For paired-end reads in bbmap operations: 'merge' or 'interleave'
"

while getopts "ho:t:" opts 
do
    case "${opts}" in
        h) echo "${usage}"; exit;;
        o) out_dir="${OPTARG}";;
        t) threads="${OPTARG}";;
        m) merge="${OPTARG}";; #NOT YET IMPLEMENTED
        [?]) echo "invalid input param"; exit 1;;
    esac
done

#################
## MAIN SCRIPT ##
#################

mkdir "${out_dir}"
mkdir "${out_dir}/merged_reads/"

# for file in "${input}/*${file_ending1}"; do
#         read1=$file
#         read2=${file%"1.fq.gz"}"2.fq.gz"

#         if test -f "$read2"; then
#             ${SCRIPT_DIR}/bbmap_pipeline.sh $read1 $read2 ${out_dir}/merged_reads/${file%"1.fq.gz"}not_merged.fastq 
#         fi     
#         rm ./tmp.*   
# done

mkdir "${out_dir}/subsampled_reads/"

for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 > 10485760' | cut -f 2 | xargs); do
   ## Initial Subampling for 28 Million Reads (Parallelized) ## 
   subsample "${file}" "${initial_sampling}" &
   pwait $(($threads-1))
done

for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 <= 10485760' | cut -f 2 | xargs); do
    ## Creating symlinks for small files ##
   ln --symbolic $(realpath "${file}") "${out_dir}/subsampled_reads/"
done

mkdir "${out_dir}/consult_out/"

${SCRIPT_DIR}/CONSULT-II/consult_search -i "${SCRIPT_DIR}/CONSULT-II/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL" -o "${out_dir}/consult_out/" \
    --query-path "${out_dir}/subsampled_reads/" \ 
    --number-matches 2 \ 
    --thread-count ${threads} \ 
    --unclassified-out "${out_dir}/consult_out/"

mkdir "${out_dir}/kraken_out/"

for file in "${out_dir}/consult_out/"; do
    ${SCRIPT_DIR}/kraken2/kraken2 --db ${SCRIPT_DIR}/kraken2/kraken_human_lib/ "${file}" --threads ${threads} --unclassified-out ${out_dir}/kraken_out/${file%".fastq"}_unclassout.fastq
done 

for file in $(ls "${out_dir}/kraken_out/"); do
    ## Runs skmer for all files
    run_skmer "${out_dir}/subsampled_reads/${file}" "${out_dir}"
done

mkdir "${out_dir}/respect_data/"

echo -e "Input\tread_length" > "${out_dir}/respect_data/hist_info.txt"
for directory in $(find "${out_dir}/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
    ## Prepares respect data directory (links skmer's .hist files)
    ## creates the hist_info.txt file that contains average read lengths
    file=${directory##*/}
    read_len=$(get_read_length "${directory}/${file}.dat")
    echo -e "${file}.hist\t${read_len}" >> "${out_dir}/respect_data/hist_info.txt"
    ln --symbolic "$(realpath ${directory}/${file}.hist)" "${out_dir}/respect_data/"
done

respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads 2
coverages=$(check_coverage ${out_dir}/estimated-parameters.txt)

for sample in $coverages; do
    genome=${sample%.hist,*}
    rm "${out_dir}/skmer_library/${genome}/${genome}.hist"
    rm "${out_dir}/skmer_library/${genome}/${genome}.dat"
    rm "${out_dir}/skmer_library/${genome}/${genome}.msh"
    rmdir "${out_dir}/skmer_library/${genome}"
    rm "${out_dir}/respect_data/${genome}.hist"
    rm "${out_dir}/subsampled_reads/${genome}.f*" 

    calc_new_subsample "${sample}" &
    pwait $(($threads-1))
done

${SCRIPT_DIR}/CONSULT-II/consult_search -i "${SCRIPT_DIR}/CONSULT-II/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL" -o "${out_dir}/consult_out/" \
    --query-path "${out_dir}/subsampled_reads/" \ 
    --number-matches 2 \ 
    --thread-count ${threads} \ 
    --unclassified-out "${out_dir}/consult_out/"

mkdir "${out_dir}/kraken_out/"

for file in "${out_dir}/consult_out/"; do
    ${SCRIPT_DIR}/kraken2/kraken2 --db ${SCRIPT_DIR}/kraken2/kraken_human_lib/ "${file}" --threads ${threads} --unclassified-out ${out_dir}/kraken_out/${file%".fastq"}_unclassout.fastq
done 

if [ "$(ls ${out_dir}/skmer_library/)" == "CONFIG" ]; then 
    rm "${out_dir}/skmer_library/CONFIG"
    rmdir "${out_dir}/skmer_library/"
fi

for file in $(ls "${out_dir}/subsampled_reads/"); do
    if ! $(ls ${out_dir}/skmer_library/ | grep -q ${file%.*}); then
        run_skmer "${out_dir}/subsampled_reads/${file}" "${out_dir}"
    fi
done

echo -e "Input\tread_length" > "${out_dir}/respect_data/hist_info.txt"
for directory in $(find "${out_dir}/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
    ## Prepares respect data directory (links skmer's .hist files)
    ## creates the hist_info.txt file that contains average read lengths
    file=${directory##*/}
    read_len=$(get_read_length "${directory}/${file}.dat")
    echo -e "${file}.hist\t${read_len}" >> "${out_dir}/respect_data/hist_info.txt"
    ln --symbolic "$(realpath ${directory}/${file}.hist)" "${out_dir}/respect_data/"
done

respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads 2
