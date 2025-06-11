import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_pav_sequence(genome_sequence, pav_coordinates, output_fasta, output_coordinates, min_length):
    
    fasta_file = SeqIO.to_dict(SeqIO.parse(genome_sequence, "fasta"))
    
    
    with open(output_coordinates, 'w') as coords, \
         open(pav_coordinates, "r") as missing_regions, \
         open(output_fasta, "w") as pav_sequences:
    
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
                
                seg_start = None  
                for i, base in enumerate(subsequence):
                    if base.upper() == 'N':
                     
                        if seg_start is not None:
                            seg_end = i - 1 
                            seg_length = seg_end - seg_start + 1
                            if seg_length >= min_length:
                                
                                out_start = start_seq + seg_start + 1  
                                out_end = start_seq + seg_end + 1
                                record_id = f"{fasta_id}__{out_start}__{out_end}"
                               
                                coords.write(f"{fasta_id}\t{out_start}\t{out_end}\n")
                             
                                segment = subsequence[seg_start:seg_end+1]
                                record = SeqRecord(Seq(segment), id=record_id, description="")
                                SeqIO.write(record, pav_sequences, "fasta")
                            seg_start = None  
                    else:
                        
                        if seg_start is None:
                            seg_start = i
                
                if seg_start is not None:
                    seg_end = len(subsequence) - 1
                    seg_length = seg_end - seg_start + 1
                    if seg_length >= min_length:
                        out_start = start_seq + seg_start + 1
                        out_end = start_seq + seg_end + 1
                        record_id = f"{fasta_id}__{out_start}__{out_end}"
                        coords.write(f"{fasta_id}\t{out_start}\t{out_end}\n")
                        segment = subsequence[seg_start:seg_end+1]
                        record = SeqRecord(Seq(segment), id=record_id, description="")
                        SeqIO.write(record, pav_sequences, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(1)
    
    genome_fasta = sys.argv[1]
    pav_coordinates = sys.argv[2]
    output_fasta = sys.argv[3]
    coordinates_out = sys.argv[4]
    threshold = int(sys.argv[5])
    
    get_pav_sequence(genome_fasta, pav_coordinates, output_fasta, coordinates_out, threshold)
