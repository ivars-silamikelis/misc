from BCBio import GFF
from Bio import SeqIO
import sys
fname=sys.argv[1]
in_file = fname
out_file = fname[:-3]+".gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
   
in_handle.close()
out_handle.close()
