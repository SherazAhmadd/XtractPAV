#-------------------------------------------------------------------------------
# Script for filtered PAV Extraction
# Authors:
#    Rana Sheraz Ahmad           : ranasheraz.202101902@gcuf.edu.pk
#    Mr. Muhammad Sadaqat        : muhammad.sadaqat@univ-rennes.fr
#    Dr. Muhammad Tahir ul Qamar : tahirulqamar@gcuf.edu.pk
# Created Time: 
#        Sat May 15 19:31:14 2025
# Version: 
#        1.0 [ initial release ]
#-------------------------------------------------------------------------------
from Bio import SeqIO
import sys
import re

def filteration(blast_out, missing_seq, qcords, unmap, coverage, identity):
    fasta_read = SeqIO.to_dict(SeqIO.parse(missing_seq, 'fasta'))

    with open(qcords, 'w') as query_coords, \
         open(qcords, 'w') as query_coords, \
         open(blast_out, 'r') as aln:

        seen = {}
        for line in aln:
            stripped = line.lstrip()
            if not stripped or stripped.startswith('#'):
                continue

            cols = re.split(r'\s+', stripped)
            if len(cols) < 10:
                continue

            qid = cols[0]
            sim = float(cols[2])
            qstart, qend = int(cols[6]), int(cols[7])
            if qid in seen:
                continue
            seen[qid] = True


            seq_len = len(fasta_read[qid].seq)
            aln_len = abs(qend - qstart) + 1
            cov = aln_len / seq_len

            if cov < float(coverage) and sim < float(identity):
                parts = qid.split('__', 2)
                if len(parts) == 3:
                    query_coords.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
                else:
                    query_coords.write(qid + "\n")

    mapped_ids = set()
    with open(blast_out, 'r') as f:
        for line in f:
            stripped = line.lstrip()
            if not stripped or stripped.startswith('#'):
                continue
            cols = re.split(r'\s+', stripped)
            if len(cols) < 1:
                continue
            mapped_ids.add(cols[0])


    unmapped_ids = set()
    for seq_id in fasta_read.keys():
        if seq_id not in mapped_ids:
            unmapped_ids.add(seq_id)

    with open(unmap, 'w') as no_hit:
        for seq_id in sorted(unmapped_ids):
            unmapped_coords = seq_id.split('__',2)
            no_hit.write(f"{unmapped_coords[0]}\t{unmapped_coords[1]}\t{unmapped_coords[2]}\n")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit("error")
    else:
        blast_input         = sys.argv[1]
        missing_sequence    = sys.argv[2]
        querycoords         = sys.argv[3]
        unmapped_coords     = sys.argv[4]
        coverage_threshold  = sys.argv[5]
        similarity_threshold= sys.argv[6]

        filteration(
            blast_input,
            missing_sequence,
            querycoords,
            unmapped_coords,
            coverage_threshold,
            similarity_threshold
        )
