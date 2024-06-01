#!/bin/bash

#set -x
#set -e

###############
## FUNCTIONS ##
###############

function subsample {
    SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    local in_file="${1}"
    local sample_size="${2}"
    local out_dir="${3}"
	
    echo "${in_file}" "${sample_size}"

    seqtk sample -s${RANDOM_SEED} -2 "${in_file}" "${sample_size}" > "${out_dir}/${in_file##*/}"
}

function calc_new_subsample {
    local sampling_data="${1}"
    local out_dir="${2}"
    local target_cov="${3}"

    local genome=$(echo $sampling_data | cut -f1 -d ',')
    local base_reads=$(echo $sampling_data | cut -f4 -d ',')
    local new_reads=$(bc <<< "${target_cov}*${base_reads}")

    subsample "${genome}" "${new_reads}" "${out_dir}"
}

function run_skmer {
    local in_file="${1}"
    local output="${2}"
    local SKMER_PROCESSORS="${3}"

    if [ -d "${output}/skmer_library" ]; then
        skmer query -a "${in_file}" "${output}/skmer_library" -p "${SKMER_PROCESSORS}" -o "${output}/temp_q"
	rm "${output}/temp_q.txt"
    else
        temp_input=$(mktemp -d "${output}/temp_in.XXXXXX")
        ln -s $(realpath "${in_file}") "${temp_input}"

        skmer reference "${temp_input}" -l "${output}/skmer_library" -p "${SKMER_PROCESSORS}" -o "${output}/temp_q"
        rm -r "${temp_input}"
	rm ${output}/temp_q*
    fi
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

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

usage="bash ${BASH_SOURCE[0]} -h [input] [-o output directory] [-t threads]

Runs nuclear read processing pipeline on a batch of merged and decontaminated reads in reference to a constructed library:
    
    Arguments:
    -h          Display this help message and exit.
    -i		    Path to INPUT directory.
    -o          Path to directory of pipeline's OUTPUT. [Default = "./OUT_fast_skims_pipeline"]
    -t          Threads to be used by most software in this pipeline (bbmap, seqtk sample, RESPECT). [Default = 2]
    -m          (T or F) Boolean, tells pipeline whether to merge or interleave paired-end reads. [Default = T] 
    -s          Size of initial sample in number of reads. [Default = 30000000]
    -c          Target coverage for subsampling. [Default = 4]
    -d          Sets top and bottom deviation thresholds for coverage (+ and - from target coverage). [Default = 1] 
    -p          Number of processes used by Skmer (large numbers of processes impacts memory). [Default = 2]
"
## TODO: Implement different number of skmer threads, post-processing pipelines, custom decontmination directories.

while getopts ":hi:o:t:c:p:f:r:s:" opts 
do
    case $opts in
        h) echo "${usage}"; exit;;
	i) INPUT="${OPTARG}" ;;
        o) OUTPUT_DIRECTORY="${OPTARG}";;
        t) NUM_THREADS="${OPTARG}";;
        c) TARGET_COV="${OPTARG}";;
        p) SKMER_PROCESSORS="${OPTARG}";;
	f) READ_1="${OPTARG}";;
	r) READ_2="${OPTARG}";;
	s) RANDOM_SEED="${OPTARG}";;
        [?]) echo "invalid input param"; exit 1;;
    esac
done

####################
## MAIN VARIABLES ##
####################

[[ -z ${INPUT} ]] && echo "Input directory not given." && exit 1
OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY:-./fast-skims_results}
NUM_THREADS=${NUM_THREADS:-1}
TARGET_COV=${TARGET_COV:-4}
SKMER_PROCESSORS=${SKMER_PROCESSORS:-1}
READ_1=${READ_1:-1.fq}
READ_2=${READ_2:-2.fq}
RANDOM_SEED=${RANDOM_SEED:-100}

#################
## MAIN SCRIPT ##
#################

mkdir "${OUTPUT_DIRECTORY}"

exec 3>&1 1>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1
set -x

###########
## SKMER ##
###########

echo "OBTAINING RAW GENOME DISTANCES..."

for file in $(ls "${INPUT}"); do
	run_skmer "${INPUT}/${file}" "${OUTPUT_DIRECTORY}" "${SKMER_PROCESSORS}"
done
skmer distance "${OUTPUT_DIRECTORY}/skmer_library/" -p "${NUM_THREADS}" -o "${OUTPUT_DIRECTORY}/distance_matrix"

rm -r "${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/"

grep "" "${OUTPUT_DIRECTORY}/skmer_library/"*/*.dat | sed -e "s/:/\t/g" -e "s/^library.//" -e "s:/[^/]*.dat\t:\t:" | sort -k2 > "${OUTPUT_DIRECTORY}/skmer_stats.tsv"

#############################
## SUBSAMPLING USING SEQTK ##
#############################

echo "SUBSAMPLING DATA..."

paste <(realpath "${INPUT}"* | parallel -j "${NUM_THREADS}" wc --lines {} | awk '{x=$1/4; printf "%s\t%.0f\n", $2, x}' | sort) \
	<(grep 'coverage' "${OUTPUT_DIRECTORY}/skmer_stats.tsv" | sort | cut  -f3) |\
	awk '{x=$2/$3; printf "%s\t%s\t%s\t%.0f\n", $1, $2, $3, x}' > "${OUTPUT_DIRECTORY}/read_counts.tsv"

export -f subsample
export -f calc_new_subsample

for coverage in ${TARGET_COV}; do
	mkdir --parents "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/"
	parallel -j "${NUM_THREADS}" calc_new_subsample {} "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/" "${coverage}" ::: "$(cat ${OUTPUT_DIRECTORY}/read_counts.tsv | tr '\t' ',')"
done

#######################################################
## PERFORMING SKIMMING OPERATIONS ON SUBSAMPLED DATA ##
#######################################################

echo "PERFORMING SKIMMING OPERATIONS ON SUBSAMPLED DATA..."

for coverage in ${TARGET_COV}; do
	for file in $(ls "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/"); do
		run_skmer "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/${file}" "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/" "${SKMER_PROCESSORS}"
	done
	skmer distance  "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/skmer_library/" -o "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/distance_matrix"

	mkdir --parents "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/"
	echo -e "Input\tread_length" > "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/hist_info.txt"
	for directory in $(find "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
    		file=${directory##*/}
    		read_len=$(get_read_length "${directory}/${file}.dat")
		
  		echo -e "${file}.hist\t${read_len}" >> "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/hist_info.txt"

    		ln --symbolic "$(realpath ${directory}/${file}.hist)" "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/"
	done


	respect --input-directories "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/" \
		--info-file "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/respect/hist_info.txt" \
		--output-directory "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/" \
		--spectra-output-size 50 \
		--iterations 1000 \
		--threads "${NUM_THREADS}"
done



