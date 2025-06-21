#-------------------------------------------------------------------------------
# Script for extraction of genic PAVs
# Authors:
#    Rana Sheraz Ahmad           : ranasheraz.202101902@gcuf.edu.pk
#    Mr. Muhammad Sadaqat        : muhammad.sadaqat@univ-rennes.fr
#    Dr. Muhammad Tahir ul Qamar : tahirulqamar@gcuf.edu.pk
# Created Time: 
#        Sat May 27 17:51:35 2025
# Version: 
#        1.0 [ initial release ]
#-------------------------------------------------------------------------------
import sys
import re

def extention(gff_file, pav_coordinates, output_coords, output_gff ):
    final_coords = open(output_coords,'w')
    final_gff = open(output_gff,'w')
    gff_data = []
    with open(gff_file,'r') as annot:
        for line in annot:
            row = re.split(r'\s+', line.strip())  
            if len(row) >= 5 and row[2] == "gene":
                gff_data.append((row[0],int(row[3]),int(row[4]),line))

    with open(pav_coordinates,'r') as coords:
        for line in coords:
            coord_sep = re.split(r'\s+', line.strip()) 
            if len(coord_sep) < 3:
                continue
            organism_name = coord_sep[0]
            start = int(coord_sep[1])
            end = int(coord_sep[2])
            found_match = False


            for gene_organism, gene_start, gene_end, gff_row in gff_data:
                if gene_organism == organism_name and gene_start <= start <= gene_end and gene_start <= end <= gene_end:
                    final_coords.write(f"{organism_name}\t{gene_start}\t{gene_end}\n")
                    final_gff.write(gff_row)
                    found_match = True
                    break

            if not found_match:
                final_coords.write(f"{organism_name}\t{start}\t{end}\n")

    final_coords.close()
    final_gff.close()
    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(1)
    else:
        gff_file = sys.argv[1]
        coordinates_file = sys.argv[2]
        output_coordinates = sys.argv[3]
        output_gff3 = sys.argv[4]
        extention(gff_file, coordinates_file, output_coordinates, output_gff3)

  
