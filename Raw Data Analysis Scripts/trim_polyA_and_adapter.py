 #split sequencing results by their indices
import sys
from analysis_helper_functions import make_save_name

file = sys.argv[1]
trim_seq = 'AAAAAAA'
new_file = sys.argv[2]
with open(file,'r') as fastq:
    with open(new_file,'w') as trimmedFile:
        #go through the fastq file, trim the trim seq, and write to output files
        l = fastq.readline()
        count = 0
        total = 0
        while l != '':
            l_1 = l
            l = fastq.readline()
            trimmed_l = l[:l.find(trim_seq)].rstrip(trim_seq)
            dif = len(l) - len(trimmed_l)
            total = total + 1
            if len(trimmed_l) > 10:
                trimmedFile.write(l_1)
                trimmedFile.write(trimmed_l + '\n')
                trimmedFile.write(fastq.readline())
                l = fastq.readline()
                trimmed_l = l[:len(l)-dif]
                trimmedFile.write(trimmed_l + '\n')
            else:
                count = count + 1
                l = fastq.readline()
                l = fastq.readline()
            l = fastq.readline()
    print("skipped {} reads because of insufficient sequence length ({} percent)".format(count, int(1000*count/total)/10.0))