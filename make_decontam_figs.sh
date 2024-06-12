#!/bin/bash

TAXONOMIC_RANK=${1:-'order'}
KRANK_REPORT_DIR=${2:-"./"}
OUT_DIR=${3:-"./"}

RESULTS='pdf(file="'${OUT_DIR}"/"${TAXONOMIC_RANK}'-decontamination_figures.pdf")'

for ABUNDANCE_PROFILE in $(realpath "$KRANK_REPORT_DIR"/abundance*); do
	RESULTS+='
	contam_data = read.csv("'${ABUNDANCE_PROFILE}'", sep = "\t")
	contam_data <- subset(contam_data, contam_data$RANK == "'${TAXONOMIC_RANK}'")
	barplot(contam_data$READ_COUNT, 
		names.arg = contam_data$TAXON_NAME, 
		col = "lightblue", 
		main = "Read Counts of Taxon Names (Rank: '${TAXONOMIC_RANK}')", 
		ylab = "log(Read Count)",
		las = 2, 
		cex.names = 0.4, 
		log = "y", 
		ylim = c(1,3*max(contam_data$READ_COUNT)))'
done 

printf "%s" "${RESULTS[@]}" | R --vanilla

