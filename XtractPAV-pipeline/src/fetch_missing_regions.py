#-------------------------------------------------------------------------------
# Script for Raw PAV Extraction
# Authors:
#    Rana Sheraz Ahmad           : ranasheraz.202101902@gcuf.edu.pk
#    Mr. Muhammad Sadaqat        : muhammad.sadaqat@univ-rennes.fr
#    Dr. Muhammad Tahir ul Qamar : tahirulqamar@gcuf.edu.pk
# Created Time: 
#        Sat May 03 05:11:47 2025
# Version: 
#        1.0 [ initial release ]
#-------------------------------------------------------------------------------

import sys

def ignoring_metadata(alignment_file):

    data = []
    with open(alignment_file, "r") as coords:
        for line in coords:
            if line.startswith("NUCMER") or line.startswith("=") or not line.strip():
                continue
            line = line.strip()
            parts = [part.strip() for part in line.split("|")]
            row = [item for part in parts for item in part.split()]
            if len(row) == 13:
                data.append(row)
    return data

def mark_covered_positions(alignment_data):

    query_coverage = {}
    query_lengths = {}

    for row in alignment_data:


        query_start, query_end = sorted([int(row[2]), int(row[3])])

        query_chr = row[12]


        query_len = int(row[8])
        

        if query_chr not in query_coverage:
            query_coverage[query_chr] = set()
            query_lengths[query_chr] = query_len

        for pos in range(query_start, query_end + 1):
            query_coverage[query_chr].add(pos)

    return  query_coverage,  query_lengths

def find_missing_regions(coverage, chrom_length, threshold):
    
    missing_regions = []
    start = None
    for pos in range(1, chrom_length + 1):
        if pos not in coverage:
            if start is None:
                start = pos  
        else:
            if start is not None:
                end = pos - 1
                if (end - start + 1) >= threshold:
                    missing_regions.append((start, end))
                start = None

    if start is not None:
        end = chrom_length
        if (end - start + 1) >= threshold:
            missing_regions.append((start, end))
    return missing_regions

def write_missing_regions(alignment_file, query_outfile, threshold):
    data = ignoring_metadata(alignment_file)
    query_cov,  query_lengths = mark_covered_positions(data)
    
    with open(query_outfile, "w") as query_f:



        for chrom in sorted(query_lengths.keys()):
            chrom_length = query_lengths[chrom]
            missing = find_missing_regions(query_cov.get(chrom, set()), chrom_length, threshold)
            for start, end in missing:
                query_f.write(f"{chrom}\t{start}\t{end}\n")
            print(f"Chromosome {chrom}: Found {len(missing)} missing regions for query.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(1)
        
    alignment_file = sys.argv[1]
    query_output = sys.argv[2]
    threshold = int(sys.argv[3])
    write_missing_regions(alignment_file, query_output, threshold)
