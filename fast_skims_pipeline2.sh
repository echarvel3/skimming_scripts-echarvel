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
    local subsampled_dir="${4}"
    local subsampled_out="subsampled_${in_file##*/}"
	
    echo "${in_file}" "${sample_size}"
    seqtk sample -s100 -2 "${in_file}" "${sample_size}" > "${out_dir}/${subsampled_dir}/${subsampled_out}"
}

function calc_new_subsample {
    local sample="${1}"
    local out_dir="${2}"
    local mean_cov="${3}"
    local initial_sampling="${4}"
    local input="${5}"
    local genome=${sample%.hist,*}
    local cov=${sample#*,}
    local optimal_samp=$(bc <<< "${mean_cov}*${initial_sampling}/${cov}")

    if echo ${genome} | grep -q "subsampled_"; then
        subsample "${input}/${genome#unclassified-seq_subsampled_*}.fastq" "${optimal_samp}" "${out_dir}" "subsampled_reads2"
    fi
}

function run_skmer {
    local in_file="${1}"
    local skmer_lib="${2}"

    if [ -d "${out_dir}/skmer_library" ]; then
        skmer query -a "${in_file}" "${skmer_lib}/skmer_library" -p "${procs}" -o "${out_dir/"temp_q"}"
    else
        temp_input=$(mktemp -d "${skmer_lib}/temp_in.XXXXXX")
        ln -s $(realpath "${in_file}") "${temp_input}"

        skmer reference "${temp_input}" -l "${skmer_lib}/skmer_library" -p "${procs}" -o "${out_dir/"temp_q"}"
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

#input=$1
def_out_dir="./OUT_fast_skims_pipeline"
def_threads=2
def_merge="T"
def_mean_cov=4
def_cov_dev=1
def_sampling=32000000
def_procs=2

file_ending1="1.fq.gz"

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

while getopts ":hi:o:t:m:s:c:d:p:" opts 
do
    case $opts in
        h) echo "${usage}"; exit;;
	    i) INPUT="${OPTARG}" ;;
        o) out_dir="${OPTARG}";;
        t) threads="${OPTARG}";;
        m) merge="${OPTARG}";;
        s) initial_sampling="${OPTARG}";;
        c) mean_cov="${OPTARG}";;
        d) cov_dev="${OPTARG}";;
        p) procs="${OPTARG}";;
        [?]) echo "invalid input param"; exit 1;;
    esac
done

# setting default values...
[[ -z ${INPUT} ]] && echo "Input directory not given." && exit 1
OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY:-./fast-skims_results}
NUM_THREADS=${NUM_THREADS:-1}
INIT_DOWNSAMPLE=${INIT_DOWNSAMPLE:-NONE}
TARGET_COV=${TARGET_COV:-4}
SKMER_PROCESSORS=${SKMER_PROCESSORS:-1}
READ_1=${READ_1:-1.fq}
READ_2=${READ_2:-2.fq}
RANDOM_SEED=${RANDOM_SEED:-100}

#################
## MAIN SCRIPT ##
#################

mkdir "${out_dir}"
mkdir --parents "${out_dir}/krank_output/krank_reports/"
mkdir --parents "${out_dir}/krank_output/decontaminated_files/"

for file in $(realpath ${input}); do
	echo -e ${file##*/}\t${file}
done > ${OUTPUT_DIRECTORY}/input_map.tsv

${SCRIPT_DIR}/KRANK/krank query \
	-l ${LIBDIR} \
	-q ${OUTPUT_DIRECTORY}/input_map.tsv \
	--max-match-distance 5 \
	--total-vote-threshold 0.03 \
	--num-threads ${NUM_THREADS} \
	--output-dir "${OUTPUT_DIRECTORY}/krank_output/krank_reports/" && echo "Decontamination step is done!"

for file in "${input}/*${file_ending1}"; do
	line_list=""
	read1=${file##*/}
	read2=${read1%${READ_1}}${READ_2}

	line_list=$(cat ${OUTPUT_DIRECTORY}/krank_output/krank_reports/classificaton-info-${read1%.*} | grep "U"  | cut -f1 | xargs -I {} grep {} -n ${input}/${read1} | cut -f1 -d ':')
	line_list+=" "$(cat ${OUTPUT_DIRECTORY}/krank_output/krank_reports/classificaton-info-${read2%.*} | grep "U"  | cut -f1 | xargs -I {} grep {} -n ${input}/${read2} | cut -f1 -d ':')
	sed_command="sed -e '"
	for line in $(echo ${line_list} | tr ' ' '\t' | uniq | xargs); do
		for  i in {0..3}; do
			sed_command+="${}d;"
		done
	done

	eval "${sed_command}' ${input}/${read1}" > "${out_dir}/krank_output/decontaminated_files/${read1}"
	eval "${sed_command}' ${input}/${read2}" > "${out_dir}/krank_output/decontaminated_files/${read2}"
done

exit 0

for line in $(cat ./classification_info-query | grep "U"  | cut -f1 | xargs -I {} grep {} -n ./query.fq | cut -f1 -d ':'); do sed_command+="${line}d;" ; done



for file in "${input}/*${file_ending1}"; do
#        read1=$file
#        read2=${file%"1.fq.gz"}"2.fq.gz"
#
#        if test -f "$read2"; then
#            ${SCRIPT_DIR}/bbmap_pipeline.sh $read1 $read2 ${out_dir}/merged_reads/${file%"1.fq.gz"}not_merged.fastq 
#        fi     
#        rm ./tmp.*   
#done

mkdir "${out_dir}/subsampled_reads/"

export -f subsample
du --summarize --block-size=1K "${input}/"* \
	| awk '$1 > 10485760' \
	| cut -f 2 \
	| parallel -j "${threads}" subsample {} "${initial_sampling}" "${out_dir}" "subsampled_reads"

for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 <= 10485760' | cut -f 2 | xargs); do
   # Creating symlinks for small files ##
 ln --symbolic $(realpath "${file}") "${out_dir}/subsampled_reads/"
done

mkdir "${out_dir}/consult_out/"

${SCRIPT_DIR}/CONSULT-II/consult_search \
  -i "${SCRIPT_DIR}/CONSULT-II/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL" \
  -o "${out_dir}/consult_out/" \
  --query-path "${out_dir}/subsampled_reads/" \
  --number-of-matches 2 \
  --thread-count ${threads} \
  --unclassified-out "${out_dir}/consult_out/"

mkdir "${out_dir}/kraken_out/"

for file in $(ls ${out_dir}/consult_out/); do
   ${SCRIPT_DIR}/kraken2/kraken2 \
	    --db ${SCRIPT_DIR}/kraken2/kraken_human_lib/ \
	    "${out_dir}/consult_out/${file}" \
	    --threads ${threads} \
	    --unclassified-out ${out_dir}/kraken_out/${file%".fastq"}.fastq
done 
 
for file in $(ls "${out_dir}/kraken_out/"); do
   ## Runs skmer for all files
   run_skmer "${out_dir}/subsampled_reads/${file}" "${out_dir}"
done

mkdir "${out_dir}/respect_data/"

echo -e "Input\tread_length" > "${out_dir}/respect_data/hist_info.txt"
for directory in $(find "${out_dir}/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
    # Prepares respect data directory (links skmer's .hist files)
    # creates the hist_info.txt file that contains average read lengths
   file=${directory##*/}
   read_len=$(get_read_length "${directory}/${file}.dat")
   echo -e "${file}.hist\t${read_len}" >> "${out_dir}/respect_data/hist_info.txt"
   ln --symbolic "$(realpath ${directory}/${file}.hist)" "${out_dir}/respect_data/"
done

#respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads "${threads}"
coverages=$(check_coverage ${out_dir}/estimated-parameters.txt)

for sample in $coverages; do
    genome=${sample%.hist,*}
    echo "genome: $genome"
    rm "${out_dir}/skmer_library/${genome}/${genome}.hist"
    rm "${out_dir}/skmer_library/${genome}/${genome}.dat"
    rm "${out_dir}/skmer_library/${genome}/${genome}.msh"
    rmdir "${out_dir}/skmer_library/${genome}"
    rm "${out_dir}/respect_data/${genome}.hist"
    rm "${out_dir}/subsampled_reads/${genome#"unclassified-seq_"}.fastq"
    rm "${out_dir}/consult_out/${genome}"
    rm "${out_dir}/kraken_out/${genome}.fastq"
done

mkdir "${out_dir}/subsampled_reads2/"

export -f calc_new_subsample
parallel -j "${threads}" calc_new_subsample {} "${out_dir}" "${mean_cov}" "${initial_sampling}" "${input}" ::: "${coverages}"

mkdir "${out_dir}/consult_out2/"
${SCRIPT_DIR}/CONSULT-II/consult_search \
   -i "${SCRIPT_DIR}/CONSULT-II/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL" \
   -o "${out_dir}/consult_out2/" \
   --query-path "${out_dir}/subsampled_reads2/" \
   --number-of-matches 2 \
   --thread-count ${threads} \
   --unclassified-out "${out_dir}/consult_out2/"

mkdir "${out_dir}/kraken_out2/"
for file in $(ls "${out_dir}/consult_out2/"); do
    ${SCRIPT_DIR}/kraken2/kraken2 \
	    --db ${SCRIPT_DIR}/kraken2/kraken_human_lib/ \
	    "${out_dir}/consult_out2/${file}" \
	    --threads ${threads} \
	    --unclassified-out ${out_dir}/kraken_out2/${file%".fastq"}.fastq
done 

if [ "$(ls ${out_dir}/skmer_library/)" == "CONFIG" ]; then 
   rm "${out_dir}/skmer_library/CONFIG"
   rmdir "${out_dir}/skmer_library/"
fi

for file in $(ls "${out_dir}/kraken_out2/"); do
   if ! $(ls ${out_dir}/skmer_library/ | grep -q ${file%.*}); then
       run_skmer "${out_dir}/kraken_out2/${file}" "${out_dir}"
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

respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads "${threads}"

