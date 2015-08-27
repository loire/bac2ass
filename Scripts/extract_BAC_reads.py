#!/usr/bin/python
import sys,re
samfile = open(sys.argv[1],'r')

# Second arguments deals with paired reads specification. Add a tag to precise if it's the forward or the reverse strand.
# The tag can be found as the second field of the first line of a fastq file, ususally (illumina) in the form of:
# Forward: 1:N:0:0 
# Reverse: 2:N:0:0

if len(sys.argv)>2:
    toadd = " "+sys.argv[2]
else:
    toadd = ""
# Creating output file name based on input
name = ".".join(sys.argv[1].split(".")[:-1])+".vector_filtered.fastq"
fout = open(name,"w")


# Regular expression used to parse cigar code from sam
cigMre = re.compile('[0-9]+M')
cigSre = re.compile('[0-9]+S')



# Parse cigar string to get length of aligment (count M) and return span of the aligned part of the read for clipping.
# Must watch for indel/ soft and hard clipping in order to be sure to keep the "good" part (the one that do NOT match). internal M must be avoided for example.  
def parse_cig(c,l):
    # c is the cigar string, l is the length of the read (necessary to get position to keep
    cigMiter = re.finditer(cigMre,c)                   
    alig_len = 0
    cpt = 0
    for match in cigMiter:
        cpt+=1
        alig_len+=int(match.group()[:-1])   # int(match.group()[:-1]) is the number before the M character in the cigar string
    # if more than one matching position on the alignment (insertion/deletion, weird stuff, we will discard the read)
    if cpt>1:
        return (0,0)
    # find soft clipped position:
    cigSiter = re.finditer(cigSre,c)
    pos_keep=[]
    for match in cigSiter:
        # if soft clipped postion is at the beginning of the read:
        if match.span()[0]==0:
            pos_keep=[0,int(match.group()[:-1])]  # int(match.group()[:-1]) is the number before the S character in the cigar string
        else:    
            #if soft clipped position is at the end of the read
            if match.span()[1]==len(c):
                pos_keep=[l-int(match.group()[:-1]),l]

    return (float(alig_len),pos_keep)

def store_reads(f,o,p):
    c = o.split()
    readname = c[0]+" "+toadd
    readseq = c[9][p[0]:p[1]]
    readqual = c[10][p[0]:p[1]]
    f.write("@"+readname+"\n")
    f.write(readseq+"\n")
    f.write("+\n")
    f.write(readqual+"\n")
    
# Treat sam line in order to choose to retain (or not) the read:

def parse_se(obj):
    c = obj.split()
#    readname = c[0]
    flag = int(c[1])
    # Test for sam flag. If> 16, means read is dodgy (mapped in several location). We will just throw it away.
    if flag > 16:
        return 0
    # else
    ref = c[2]
    readlen = len(c[9])
#    MapQ = int(c[4])
#    cig = c[5]
    # If the reads didn't map on the vector, then we keep it as such. I'm not sure that there's no more vector sequence in this read (if the vector sequence is shorter than baw seed), so a trimming of the edge will be necessary
    if ref=="*":
        store_reads(fout,obj,[0,readlen])
    # Else, the read map to the vector. Instead of throwing it away, we will examine it in order to keep the part of the read that don't map to the vector, based on the cigar sequence.
    else:
        readname = c[0]
        cig = c[5]
        # get the aligned part.
        alig_len,Spos = parse_cig(cig,readlen)
        # If dodgy read (align in more than one part, ie cigar M is splitted on the ref, disregard read (only a few are concerned)
        if alig_len == 0:
            return 0
        
        #print Spos
        #NM flag is in the 11th position
        NM = c[11]
        #edit distance: last field of NM
        eNM = int(NM.split(":")[-1])
        
        # pc_id: identity percentage of the mapped region / cov : proportion of the reads that map to the vector
        pc_id = 1-(eNM/alig_len)
        cov = alig_len/float(readlen)
        #print obj[:-1]
        #print readname,ref,alig_len,eNM,pc_id,readlen,alig_len,cov
        
        # So if the pc_id is high and at least 10 percent of the read is recoverable, 
        if pc_id >= 0.95 and cov <= 0.9:
            store_reads(fout,obj,Spos)



line = "@SQ"
while line[0:3]=="@SQ" or line[0:3]=="@PG":
    line = samfile.readline()
cpt = 0

# here for multithreading if you feel like implmenting it.
# This prog with bwa mem is a good alternative to deconseq, think about it. 
for line in samfile:
#    print line.rstrip()
    parse_se(line)
    cpt+=1
#    if cpt > 2000:
#        break
fout.close()
samfile.close()


