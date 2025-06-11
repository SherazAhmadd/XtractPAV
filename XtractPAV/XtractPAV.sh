#!/bin/bash
# Author: Rana Sheraz Ahmad [ ranasheraz.202101902@gcuf.edu.pk ]
# XtractPAV: A Presence/Absence Variation Scanning and Identification tool
# Version: 1.0

version="XtractPAV v1.0"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)
            echo -e "\n $version"
            echo ""
            echo "[  Usage  ]"
            echo ""
            echo " XtractPAV.sh --rf <reference_genome.fa> --ra <reference_genome.gff3> --qf <query_genome1.fa,query_genome2.fa> --qa <query_genome1.gff3,query_genome2.gff3> --cov <float> --sim <float> --len <int> --thr <int>"
            echo ""
            echo "[  Options  ]"
            echo ""
            echo "      --rf              Reference genome FASTA file"
            echo "      --ra              Reference genome GFF3 annotation file"
            echo "      --qf              Comma-separated query genome FASTA files"
            echo "      --qa              Comma-separated query genome GFF3 files"
            echo "      --cov             Coverage threshold (e.g., 0.9)"
            echo "      --sim             Similarity threshold (e.g., 95.0)"
            echo "      --len             Minimum length of PAVs (e.g., 100)"
            echo "      --thr             Number of threads to use"
            echo "      --help            Show this help message and exit"
            echo "      --version         Show script version and exit"
            echo ""
            echo "[  Dependencies  ]"
            echo ""
            echo "      - MUMmer4"
            echo "      - BLAST+"
            echo "      - Bedtools v2.27.1"
            echo "      - Python3"
            echo "      - Biopython==1.78"
            echo "      - Plotly==6.0.1"
            echo ""
            exit 0
            ;;
        --version)
            echo "$version"
            exit 0
            ;;
        --rf) ref_genome="$2"; shift 2 ;;
        --ra) ref_gff3="$2"; shift 2 ;;
        --qf) IFS=',' read -r -a query_genomes <<< "$2"; shift 2 ;;
        --qa) IFS=',' read -r -a query_gffs <<< "$2"; shift 2 ;;
        --cov) cov="$2"; shift 2 ;;
        --sim) sim="$2"; shift 2 ;;
        --len) th="$2"; shift 2 ;;
        --thr) thr="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
done




usage() {
    echo -e "\n$version"
    echo ""
    echo "[  Usage  ]"
    echo ""
    echo "XtractPAV.sh --rf <reference_genome.fa> --ra <reference_genome.gff3> --qf <query_genome1.fa,query_genome2.fa> --qa <query_genome1.gff3,query_genome2.gff3> --cov <float> --sim <float> --len <int> --thr <int>"
    echo ""
    echo "Run 'XtractPAV.sh --help' for more details."
    exit 1
}

if [[ -z "$ref_genome" || -z "$ref_gff3" || -z "$query_genomes" || -z "$query_gffs" || -z "$cov" || -z "$sim" || -z "$th" || -z "$thr" ]]; then
    echo "[  Error  ] Missing required arguments."
    usage
fi

clear

log_time() {
    local step_label="$1"
    local start_time=$2
    local end_time=$(date +%s)
    local elapsed_time=$((end_time - start_time))
    local elapsed_minutes=$(awk "BEGIN {printf \"%.2f\", $elapsed_time/60}")
    echo -e "\033[1;34m$step_label Execution Time: $elapsed_time seconds ($elapsed_minutes minutes)\033[0m"
}

overall_start=$(date +%s)

echo -e "\033[1;30m________________________________________________________________________________________________________________________________________________________________________\033[0m\n"
echo -e "\033[1;32m\t\t\t\t\tWelcome to the Presence/Absence Variant (PAV) Scanning Pipeline for Eukaryotes and Prokaryotes\033[0m\n"
echo -e "\033[1;30m________________________________________________________________________________________________________________________________________________________________________\033[0m\n"
echo -e "\033[1;34m----------------------------------------------------------------------------- Parameter Explanation --------------------------------------------------------------------\033[0m"
echo -e "\033[1;30m________________________________________________________________________________________________________________________________________________________________________\033[0m\n"
echo -e "\n"

echo -e "\033[1;33m\tCoverage Threshold:\033[0m This threshold is used to filter Presence/Absence Variants (PAVs) based on coverage."
echo -e "\tThe value should be between 0 and 1. \033[1;32mRecommended: 0.9\033[0m"
echo -e "\n"

echo -e "\033[1;33m\tSimilarity Threshold:\033[0m The similarity value filters the PAVs based on similarity percentage."
echo -e "\tThe value should be between 0 and 1. \033[1;32mRecommended: 95.0%\033[0m"
echo -e "\n"

echo -e "\033[1;33m\tPAV Length:\033[0m Defines the minimum base pair length for Presence/Absence Regions in the genome."
echo -e "\tThis determines how long the PAVs should be for inclusion in the analysis."
echo -e "\tThe length should be at least \033[1;32m100bp\033[0m."
echo -e "\n"

echo -e "\033[1;30m_______________________________________________________________________________________________________________________________________________________________________\033[0m\n"
echo -e "\n\033[1;34m----------------------------------------------------------------------------- Configuration Summary --------------------------------------------------------------------\033[0m\n"

echo -e "\033[1;30m_______________________________________________________________________________________________________________________________________________________________________\033[0m\n"
echo -e "\033[1;32mConfiguration Parameters\033[0m\n"
echo -e "\033[1;33m\t- Coverage Threshold fixed at \033[0m$cov"
echo -e "\033[1;33m\t- Similarity Threshold fixed at \033[0m$sim"
echo -e "\033[1;33m\t- Minimum PAV Length fixed at \033[0m$th"
echo -e "\n\033[1;31mProceeding with the pipeline execution...\033[0m"
echo -e "\033[1;30m______________________________________________________________________________________________________________________________________________________________________\033[0m"
echo -e "\n"



script_dir="$(dirname "$(realpath "$0")")"


if [[ -d "$script_dir/src" ]]; then
    main_dir="$script_dir"
else
    main_dir="$(dirname "$script_dir")"
fi


if [[ ! -d "$main_dir/src" ]]; then
    echo "Error: 'src' directory not found relative to script ($script_dir)"
    exit 1
fi


for i in "${!query_genomes[@]}"; do
    query_genome="${query_genomes[$i]}"
    query_gff3="${query_gffs[$i]}"

    query_name=$(basename "$query_genome")
    query_base_name="${query_name%.*}"
    query_output_dir="./Results/$query_base_name"


    mkdir -p "$query_output_dir"
    temp_dir="$main_dir/temp"
    mkdir -p "$temp_dir"

    echo -e "\033[1;32mStep - 1a : \033[0m Running nucmer for genome alignment for $query_base_name"
    step_start=$(date +%s)
    nucmer --prefix="$temp_dir/nucmer_output" -t "$thr" "$ref_genome" "$query_genome"
    if [[ $? -ne 0 ]]; then
        echo "nucmer execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "nucmer" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 1b : \033[0m Filtering nucmer alignment to retain best 1-to-1 matches"
    
    step_start=$(date +%s)
    delta-filter -1 "$temp_dir/nucmer_output.delta" > "$temp_dir/nucmer_output.filtered.delta"
    if [[ $? -ne 0 ]]; then
        echo "delta-filter execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "delta-filter" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 1c : \033[0m Generating alignment coordinates with show-coords "
    step_start=$(date +%s)
    show-coords -rcl "$temp_dir/nucmer_output.filtered.delta" > "$temp_dir/testing_alignment.coords"
    if [[ $? -ne 0 ]]; then
        echo "show-coords execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "show-coords" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 2 : \033[0m Extracting the missing region corresponding to the reference genome"
    step_start=$(date +%s)
    python3 $main_dir/src/fetch_missing_regions.py "$temp_dir/testing_alignment.coords" "$temp_dir/query_out.txt" $th
    if [[ $? -ne 0 ]]; then
        echo "fetch_missing_regions.py execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "fetch_missing_regions.py" $step_start
    echo -e "\n"


    echo -e "\033[1;32mStep - 3 : \033[0m Retrieving the Sequence of missing Regions from genome fasta file"
    step_start=$(date +%s)
    python3 $main_dir/src/fetch_sequence.py "$query_genome" "$temp_dir/query_out.txt" "$temp_dir/missing_seq.fa" "$temp_dir/pav_coordinates.txt" $th
    if [[ $? -ne 0 ]]; then
        echo "fetch_sequence.py execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "fetch_sequence.py" $step_start
    echo -e "\n"
    echo -n "Total number of PAVs extracted from $query_base_name: "
    # echo -e "\033[1;31m$(grep -c '>' $temp_dir/filtered_seq.fa)\033[0m"
    grep -c '>' $temp_dir/missing_seq.fa
    echo -e "\n"
    

    echo -e "Running BLAST to align missing sequences"
    step_start=$(date +%s)
    makeblastdb -in "$ref_genome" -dbtype nucl -out "$temp_dir/short_testing_db"
    if [[ $? -ne 0 ]]; then
        echo "makeblastdb execution failed for $query_base_name. Exiting."
        exit 1
    fi

    blastn -query "$temp_dir/missing_seq.fa" -db "$temp_dir/short_testing_db" -outfmt 6 -out "$temp_dir/blast_testing.txt" -evalue 1e-5 -num_threads "$thr"
    if [[ $? -ne 0 ]]; then
        echo "blastn execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "BLAST" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 4 : \033[0m Filteration process: Applying the Coverage and Similarity Criteria to Shortlist The Authentic PAVs"
    step_start=$(date +%s)
    python3 $main_dir/src/fetch_filtered_PAVs.py "$temp_dir/blast_testing.txt" "$temp_dir/missing_seq.fa" "$temp_dir/criteria_defined_PAVs.bed" "$temp_dir/No-Hit_PAVs.bed" $cov $sim
    if [[ $? -ne 0 ]]; then
        echo "fetch_filtered_PAVs.py execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "fetch_filtered_PAVs.py" $step_start

    cat "$temp_dir/criteria_defined_PAVs.bed" "$temp_dir/No-Hit_PAVs.bed" | sortBed > $temp_dir/total_pavs.bed
    echo -e "\n"
    echo -n "The  number of high-confidence PAVs with defined stringent filtering criteria in $query_base_name: "
    cat $temp_dir/criteria_defined_PAVs.bed | wc -l
    echo -n "The  number of PAVs with No Hit / unique in $query_base_name: "
    cat $temp_dir/No-Hit_PAVs.bed  | wc -l
    echo -e "\n"
    echo -n "Total number of PAVs identified after filteration: "
    cat $temp_dir/total_pavs.bed | wc -l
    echo -e "\n"

    echo -e "Extracting The Gene Annotations From The Query GFF3 File"
    step_start=$(date +%s)
    grep -P "\tgene\t" "$query_gff3" > "$temp_dir/${query_base_name}_loci.gff3"
    if [[ $? -ne 0 ]]; then
        echo "Gene annotation filtering failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "Gene annotation extraction" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 5 : \033[0m Extending The PAVs and Extracting The Genic PAVs and Their Annotation File"
    step_start=$(date +%s)
    python3 $main_dir/src/fetch_genic_PAVs.py "$temp_dir/${query_base_name}_loci.gff3" "$temp_dir/total_pavs.bed" "$temp_dir/Pav_final_coords.txt" "$query_output_dir/${query_base_name}_Genic_PAVs.gff3"
    if [[ $? -ne 0 ]]; then
        echo "fetch_genic_PAVs.py execution failed for $query_base_name. Exiting."
        exit 1
    fi
    log_time "fetch_genic_PAVs.py" $step_start
    echo -e "\n"
    # echo -n "Number of genic-PAVs found in $query_base_name: "
    # echo -e "\033[1;31m$(grep -cP "\tgene\t" $query_output_dir/${query_base_name}_Genic_PAVs.gff3)\033[0m"
    # grep -cP "\tgene\t" $query_output_dir/${query_base_name}_Genic_PAVs.gff3
    # echo -e "\n"

    echo -e "\033[1;32mStep - 6 : \033[0m Retrieving The Extended PAV Sequence To Cover The Gene From Genome File"
    step_start=$(date +%s)
    python3 $main_dir/src/fetch_genic_sequence.py "$query_genome" "$temp_dir/Pav_final_coords.txt" "$query_output_dir/${query_base_name}.fa" "$temp_dir/Coord.txt" 
    if [[ $? -ne 0 ]]; then
        echo "fetch_genic_sequence.py execution failed. Exiting."
        exit 1
    fi
    log_time "fetch_genic_sequence.py" $step_start
    echo -e "\n"

    echo -e "\033[1;32mStep - 7 : \033[0m Generating The Final Summary Report of PAV analysis"
    step_start=$(date +%s)    
    python3 $main_dir/src/fetch_summary.py "$temp_dir/total_pavs.bed" "$query_output_dir/${query_base_name}_Genic_PAVs.gff3" "$query_output_dir/${query_base_name}.html"
    if [[ $? -ne 0 ]]; then
        echo "fetch_summary.py execution failed. Exiting."
        exit 1
    fi
    log_time "fetch_summary.py" $step_start
    echo -e "\n"
    echo -e "\033[1;30m______________________________________________________________________________________________________________________________________________________________________\033[0m"
    echo -e "\n"

done


overall_end=$(date +%s)
total_time=$((overall_end - overall_start))
total_minutes=$(awk "BEGIN {printf \"%.2f\", $total_time/60}")
total_hours=$(awk "BEGIN {printf \"%.4f\", $total_time/3600}")

echo -e "\033[1;32mXtractPAV Pipeline has successfully completed the Presence/Absence Variation analysis for all provided genomes. Results, including annotated genic PAVs and the summary report, have been systematically archived in the directory: $query_output_dir\033[0m\n"
echo -e "\033[1;31mTotal Execution Time: $total_time seconds ($total_minutes minutes, $total_hours CPU hours)\033[0m"
