#!/bin/bash

#set -x
#set -e

###############
## FUNCTIONS ##
###############

function subsample {
    local in_file="${1}"
    local sample_size="${2}"
    local out_dir="${3}"
    local SCRIPT_DIR="${4}"
    
    echo "${in_file}" "${sample_size}"

    in_file1=$(mktemp "${out_dir}/tmp.XXXXXX")
    in_file2=$(mktemp "${out_dir}/tmp.XXXXXX")

    ${SCRIPT_DIR}/bbmap/reformat.sh in="${in_file}" out1="${in_file1}" out2="${in_file2}" overwrite=true

    seqtk sample -s${RANDOM_SEED} -2 "${in_file1}" "$(echo ${sample_size}/2 | bc)" > "${out_dir}/temp2-${in_file##*/}"
    seqtk sample -s${RANDOM_SEED} -2 "${in_file2}" "$(echo ${sample_size}/2 | bc)" > "${out_dir}/temp1-${in_file##*/}"
    rm "${in_file1}" "${in_file2}"
    
    ${SCRIPT_DIR}/bbmap/reformat.sh in1="${out_dir}/temp2-${in_file##*/}" in2="${out_dir}/temp1-${in_file##*/}" out="${out_dir}/${in_file##*/}" overwrite=true
    rm "${out_dir}/temp1-${in_file##*/}" "${out_dir}/temp2-${in_file##*/}"
}

function calc_new_subsample {
    local sampling_data="${1}"
    local out_dir="${2}"
    local target_cov="${3}"
    local SCRIPT_DIR="${4}"
    
    local genome=$(echo $sampling_data | cut -f1 -d ',')
    local base_reads=$(echo $sampling_data | cut -f4 -d ',')
    local new_reads=$(bc <<< "${target_cov}*${base_reads}")

    subsample "${genome}" "${new_reads}" "${out_dir}" "${SCRIPT_DIR}"
}

function run_skmer {
    local in_file="${1}"
    local output="${2}"

    if [ -d "${output}/skmer_library" ]; then
        skmer query -a "${in_file}" "${output}/skmer_library" -p "${SKMER_PROCESSORS}" -o "${output}/temp_q"
	rm ${output}/temp_q*
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
    -i	        Path to INPUT directory.
    -o          Path to directory of pipeline's OUTPUT. [Default = "./fast-skims_results"]
    -t          Threads to be used by most software in this pipeline (bbmap, seqtk sample, RESPECT). [Default = 1]
    -c          Target coverage for subsampling. [Default = 4]		
    -p          Number of processes used by Skmer (large numbers of processes impacts memory). [Default = 2]
    -f		ending format for Read 1 of paired end reads (e.g. "_R1.fastq" or ".01.fq") [Default = 1.fq]
    -r		ending format for Read 2 of paired end reads (e.g. "_R2.fastq" or ".02.fq") [Default = 2.fq]
    -s		random seed used for all software in the pipeline.
    -d		downsampling method (SKMER or RESPECT) [Default = SKMER]
"
## TODO: Implement post-processing pipelines, custom decontmination directories.

while getopts ":hi:o:t:c:p:f:r:s:d:l:" opts 
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
	d) DOWNSAMPLING_VERSION="${OPTARG}";;
	l) LIBRARIES="${OPTARG}";;
        [?]) echo "invalid input param"; exit 1;;
    esac
done

####################
## MAIN VARIABLES ##
####################

[[ -z ${INPUT} ]] && echo "Input directory not given." && exit 1
OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY:-./fast-skims_results}
mkdir "${OUTPUT_DIRECTORY}"

exec 3>&1 1>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

NUM_THREADS=${NUM_THREADS:-1}
TARGET_COV=${TARGET_COV:-4}
SKMER_PROCESSORS=${SKMER_PROCESSORS:-1}
READ_1=${READ_1:-1.fq}
READ_2=${READ_2:-2.fq}
RANDOM_SEED=${RANDOM_SEED:-100}
DOWNSAMPLING_VERSION=${DOWNSAMPLING_VERSION:-SKMER}
LIBRARIES=${LIBRARIES:-"${SCRIPT_DIR}/KRANK/lib_reps_adpt-k29_w35_h13_b16_s8 ${SCRIPT_DIR}/KRANK/pangenome-05-2024-lib_rand_free-k29_w34_h13_b16_s8"}

#################
## MAIN SCRIPT ##
#################
set -x

mkdir --parents "${OUTPUT_DIRECTORY}"

#####################
## DECONTAMINATION ##
#####################
exec 1>&3
echo "DECONTAMINATING DATA..."
exec 3>&1 
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

###########
## SKMER ##
###########

echo "OBTAINING GENOME DISTANCES FROM SEQUENCES..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

for file in $(ls "${INPUT}/"); do
	run_skmer "${INPUT}/${file}" "${OUTPUT_DIRECTORY}"
done

skmer distance "${OUTPUT_DIRECTORY}/skmer_library/" -p "${NUM_THREADS}" -o "${OUTPUT_DIRECTORY}/distance_matrix"



