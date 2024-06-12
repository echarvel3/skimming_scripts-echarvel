#!/bin/bash

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
        skmer query -a "${in_file}" "${skmer_lib}/skmer_library" -p "2"
    else
        temp_input=$(mktemp -d "${skmer_lib}/temp_in.XXXXXX")
        ln -s $(realpath "${in_file}") "${temp_input}"

        skmer reference "${temp_input}" -l "${skmer_lib}/skmer_library" -p "2"
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
#input=$1
def_out_dir="./OUT_subsample_and_estimate"
def_threads=2
def_mean_cov=4
def_cov_dev=1
def_sampling=30000000

usage="bash ${BASH_SOURCE[0]} -h [ -i input ] [ -o output directory ] [ -t threads ] [ -m merge ] [ -s initial read sampling ] [ -c target coverage ] [ -d coverage deviation ]

Runs nuclear read processing pipeline on a batch of merged and decontaminated reads in reference to a constructed library:
    
    Arguments:
    -h          Display this help message and exit.
    -i          Path to INPUT directory.
    -o          Path to directory of pipeline's OUTPUT. [Default = "./OUT_fast_skims_pipeline"]
    -t          Threads to be used by all software in this pipeline (seqtk sample, Skmer, RESPECT). [Default = 2]
    -s          Size of initial sample in number of reads. [Default = 30000000]
    -c          Target coverage for subsampling. [Default = 4]
    -d          Sets top and bottom deviation thresholds for coverage (+ and - from target coverage). [Default = 1] 
"
## TODO: Implement different number of skmer threads, post-processing pipelines, custom decontmination directories.

while getopts ":hi:o:t:m:s:c:d:" opts 
do
    case $opts in
        h) echo "${usage}"; exit;;
        i) input="${OPTARG}" ;;
        o) out_dir="${OPTARG}";;
        t) threads="${OPTARG}";;
        s) initial_sampling="${OPTARG}";;
        c) mean_cov="${OPTARG}";;
        d) cov_dev="${OPTARG}";;
        [?]) echo "invalid input param"; exit 1;;
    esac
done

# setting default values...
[[ -z $input ]] && echo "NO INPUT GIVEN"; exit 1
[[ -z $out_dir ]] && out_dir="${def_out_dir}"
[[ -z $threads ]] && threads="${def_threads}"
[[ -z $initial_sampling ]] && initial_sampling="${def_sampling}"
[[ -z $mean_cov ]] && mean_cov="${def_mean_cov}"
[[ -z $cov_dev ]] && cov_dev="${def_cov_dev}"
high_cov="$(bc <<< $mean_cov + $cov_dev)"
low_cov="$(bc <<< $mean_cov - $cov_dev)"

#################
## MAIN SCRIPT ##
#################

mkdir "${out_dir}"
mkdir "${out_dir}/subsampled_reads/"

export -f subsample
du --summarize --block-size=1K "${input}/"* \
	| awk '$1 > 10485760' \
	| cut -f 2 \
	| parallel -j "${threads}" subsample {} "${initial_sampling}" "${out_dir}" "subsampled_reads"

for file in $(du --summarize --block-size=1K "${input}/"* | awk '$1 <= 10485760' | cut -f 2 | xargs); do
    ## Creating symlinks for small files ##
   ln --symbolic $(realpath "${file}") "${out_dir}/subsampled_reads/"
done

for file in $(ls "${out_dir}/subsampled_reads/"); do
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

respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads "${threads}"
coverages=$(check_coverage ${out_dir}/estimated-parameters.txt)

for sample in $coverages; do
    genome=${sample%.hist,*}
    rm "${out_dir}/skmer_library/${genome}/${genome}.hist"
    rm "${out_dir}/skmer_library/${genome}/${genome}.dat"
    rm "${out_dir}/skmer_library/${genome}/${genome}.msh"
    rmdir "${out_dir}/skmer_library/${genome}"
    rm "${out_dir}/respect_data/${genome}.hist"
    rm "${out_dir}/subsampled_reads/${genome}.f*" 
done

export -f calc_new_subsample
parallel -j "${threads}" calc_new_subsample {} "${out_dir}" "${mean_cov}" "${initial_sampling}" "${input}" ::: "${coverages}"

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

respect -d "${out_dir}/respect_data/" -I "${out_dir}/respect_data/hist_info.txt" -o "${out_dir}" -N 10 --threads "${threads}"
