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

while getopts ":hi:o:t:c:p:f:r:s:d:" opts 
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

#################
## MAIN SCRIPT ##
#################
set -x

mkdir --parents "${OUTPUT_DIRECTORY}/krank_output/krank_reports/"
mkdir --parents "${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/"

#####################
## DECONTAMINATION ##
#####################
exec 1>&3
echo "DECONTAMINATING DATA..."
exec 3>&1 
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1


for file in $(realpath ${INPUT}/*); do
	# ...creates input map for KRANK input...
	echo -e "${file##*/}\t${file}"
done > ${OUTPUT_DIRECTORY}/input_map.tsv

# ...running KRANK decontamination...
${SCRIPT_DIR}/KRANK/krank query \
	--library-dir ${SCRIPT_DIR}/KRANK/lib_reps_adpt-k29_w35_h13_b16_s8 ${SCRIPT_DIR}/KRANK/pangenome-05-2024-lib_rand_free-k29_w34_h13_b16_s8 \
	--query-file ${OUTPUT_DIRECTORY}/input_map.tsv \
	--max-match-distance 5 \
	--total-vote-threshold 0.03 \
	--num-threads ${NUM_THREADS} \
	--output-dir "${OUTPUT_DIRECTORY}/krank_output/krank_reports/" && echo "-- KRANK run was successful!"

tmp_dir=$(mktemp --directory "${OUTPUT_DIRECTORY}/krank_output/tmp.XXXXXXX")

# ...creates decontaminated read files from KRANK output...
for file in $(realpath "${INPUT}/*${READ_1}"); do
	read1=${file##*/}
	read2=${read1%${READ_1}}${READ_2}

	cat ${OUTPUT_DIRECTORY}/krank_output/krank_reports/classification_info-${read1} | grep --perl-regexp '\tU\t'  | cut -f1 | sed -e 's/@//g' > "${tmp_dir}/unclass_readnames1.txt"
	cat ${OUTPUT_DIRECTORY}/krank_output/krank_reports/classification_info-${read2} | grep --perl-regexp '\tU\t'  | cut -f1 | sed -e 's/@//g' > "${tmp_dir}/unclass_readnames2.txt"

	${SCRIPT_DIR}/bbmap/filterbyname.sh in1=${INPUT}/${read1} in2=${INPUT}/${read2} out1="${tmp_dir}/${read1}" out2="${tmp_dir}/${read2}" threads=${NUM_THREADS} names="${tmp_dir}/unclass_readnames1.txt" include=true overwrite=true
	${SCRIPT_DIR}/bbmap/filterbyname.sh in1=${tmp_dir}/${read2} in2=${tmp_dir}/${read1} out1="${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/${read2}" out2="${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/${read1}" threads=${NUM_THREADS} names="${tmp_dir}/unclass_readnames2.txt" include=true overwrite=true

	rm ${tmp_dir}/${read2} ${tmp_dir}/${read1}
done
rm -r "${tmp_dir}"

##################################
## BBTOOLS: TRIMMING AND DEDUPE ##
##################################

exec 1>&3
echo "PERFORMING BBMAP OPERATIONS..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

mkdir "${OUTPUT_DIRECTORY}/bbmap_reads/"

for file in $(realpath "${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/*${READ_1}"); do
        read1=$file
        read2=${read1%${READ_1}}${READ_2}
	out_read=${file##*/}

        if test -f "$read2"; then
            ${SCRIPT_DIR}/bbmap_pipeline.sh "$read1" "$read2" "${OUTPUT_DIRECTORY}/bbmap_reads/${out_read%${READ_1}}.fastq" 
        fi     
done

###########
## SKMER ##
###########

exec 1>&3
echo "OBTAINING GENOME DISTANCES FROM FULL COVERAGE SEQUENCES..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

for file in $(ls "${OUTPUT_DIRECTORY}/bbmap_reads/"); do
	run_skmer "${OUTPUT_DIRECTORY}/bbmap_reads/${file}" "${OUTPUT_DIRECTORY}"
done
skmer distance "${OUTPUT_DIRECTORY}/skmer_library/" -p "${NUM_THREADS}" -o "${OUTPUT_DIRECTORY}/distance_matrix"

rm -r "${OUTPUT_DIRECTORY}/krank_output/decontaminated_files/"

grep "" "${OUTPUT_DIRECTORY}/skmer_library/"*/*.dat | sed -e "s/:/\t/g" -e "s/^library.//" -e "s:/[^/]*.dat\t:\t:" | sort -k2 > "${OUTPUT_DIRECTORY}/skmer_stats.tsv"

#############################
## RUNNING RESPECT ON DATA ##
#############################

exec 1>&3
echo "RUNNING RESPECT ON FULL COVERAGE SEQUENCES..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

mkdir --parents "${OUTPUT_DIRECTORY}/respect/"
echo -e "Input\tread_length" > "${OUTPUT_DIRECTORY}/respect/hist_info.txt"
for directory in $(find "${OUTPUT_DIRECTORY}/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
	file=${directory##*/}
	read_len=$(get_read_length "${directory}/${file}.dat")
		
	echo -e "${file}.hist\t${read_len}" >> "${OUTPUT_DIRECTORY}/respect/hist_info.txt"

	ln --symbolic "$(realpath ${directory}/${file}.hist)" "${OUTPUT_DIRECTORY}/respect/"
done

respect --input-directories "${OUTPUT_DIRECTORY}/respect/" \
	--info-file "${OUTPUT_DIRECTORY}/respect/hist_info.txt" \
	--output-directory "${OUTPUT_DIRECTORY}/" \
	--spectra-output-size 50 \
	--iterations 1000 \
	--threads "${NUM_THREADS}"

#############################
## SUBSAMPLING USING SEQTK ##
#############################

exec 1>&3
echo "SUBSAMPLING DATA..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

tmp_coverages=$(mktemp "${OUTPUT_DIRECTORY}/coverages.XXXX")
grep 'coverage' "${OUTPUT_DIRECTORY}/skmer_stats.tsv" | sort | cut  -f3 | sed -e 's/NA/1.0/g' > ${tmp_coverages}

if [[ "${DOWNSAMPLING_VERSION}" == "RESPECT" ]]; then
	cut -f1,4 estimated-parameters.txt | tail -n+2 | sort|  sed -e 's/hist\t/fastq\t/g' | cut -f2 | sed -e 's/NA/1.0/g' > ${tmp_coverages}
fi

paste <(realpath "${OUTPUT_DIRECTORY}/bbmap_reads/"* | parallel -j "${NUM_THREADS}" wc --lines {} | awk '{x=$1/4; printf "%s\t%.0f\n", $2, x}' | sort) \
	<(cat ${tmp_coverages}) |\
	awk '{x=$2/$3; printf "%s\t%s\t%s\t%.0f\n", $1, $2, $3, x}' > "${OUTPUT_DIRECTORY}/read_counts.tsv"
rm "${tmp_coverages}"

export -f subsample
export -f calc_new_subsample

for coverage in ${TARGET_COV}; do
	mkdir --parents "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/"
	parallel -j "${NUM_THREADS}" calc_new_subsample {} "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/" "${coverage}" "${SCRIPT_DIR}" ::: "$(cat ${OUTPUT_DIRECTORY}/read_counts.tsv | tr '\t' ',')"
done

#######################################################
## PERFORMING SKIMMING OPERATIONS ON SUBSAMPLED DATA ##
#######################################################

exec 1>&3
echo "PERFORMING SKIMMING OPERATIONS ON SUBSAMPLED DATA..."
exec 3>&1
exec 3>&1 1>>"${OUTPUT_DIRECTORY}/skimming-scripts.log" 2>&1

for coverage in ${TARGET_COV}; do
	for file in $(ls "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/"); do
		run_skmer "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/reads/${file}" "${OUTPUT_DIRECTORY}/subsampled_data/${coverage}x_data/"
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



