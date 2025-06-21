#-------------------------------------------------------------------------------
# Script for retrieval of genic PAVs sequence without N's
# Authors:
#    Rana Sheraz Ahmad           : ranasheraz.202101902@gcuf.edu.pk
#    Mr. Muhammad Sadaqat        : muhammad.sadaqat@univ-rennes.fr
#    Dr. Muhammad Tahir ul Qamar : tahirulqamar@gcuf.edu.pk
# Created Time: 
#        Sat May 31 01:11:39 2025
# Version: 
#        1.0 [ initial release ]
#-------------------------------------------------------------------------------
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_pav_sequence(genome_sequence, pav_coordinates, output_fasta, output_coordinates):

    fasta_file = SeqIO.to_dict(SeqIO.parse(genome_sequence, "fasta"))

    coords = open(output_coordinates,'w')

    with open(pav_coordinates, "r") as missing_regions, open(output_fasta, "w") as pav_sequences:
       for interval in missing_regions:
        interval = interval.strip()
        if not interval:
           continue
        pav = interval.split("\t")
        fasta_id = pav[0]
        start_seq = int(pav[1]) - 1
        end_seq = int(pav[2])

        if fasta_id in fasta_file:
           full_sequence = str(fasta_file[fasta_id].seq)
           subsequence = full_sequence[start_seq:end_seq]

           record_id = f"{fasta_id}__{start_seq+1}__{end_seq}"
           coords.write(f"{fasta_id}\t{start_seq+1}\t{end_seq}\n")
           record = SeqRecord(Seq(subsequence),id=record_id, description="")
        
           SeqIO.write(record, pav_sequences, "fasta")
       
if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(1)

    genome_fasta = sys.argv[1]
    pav_coordinates = sys.argv[2]
    output_fasta = sys.argv[3]
    coordinates_out = sys.argv[4]


    get_pav_sequence(genome_fasta, pav_coordinates, output_fasta, coordinates_out)
    

