import sys


def main():

    infile  = sys.argv[1]
    outfile = sys.argv[2]

    with open(infile, 'r') as f:
        with open(outfile,'w') as out:

            out.write( "track type=wiggle_0" )
            current_chrom = None
            total = 0
            reads = {}

            for line in f:

                chrom,pos,count = line.strip().split()
                pos,count = int(pos),float(count)

                # If this is the first entry I've ever read, store the chrom name
                if current_chrom is None:
                    current_chrom = chrom

                # If the current line has a different chromosome,
                # write out a new header and update the current_chrom variable
                if chrom != current_chrom:

                    if total > 0:
                        print( current_chrom,total)
                    out.write( "\n" + "variableStep chrom="+current_chrom )
                    for my_pos,my_count in sorted(reads.items(),
                                            key=lambda x: x[0]):
                        out.write("\n"+str(my_pos)+"\t"+str(my_count))

                    # Reset the variables for the new chromosome
                    current_chrom = chrom
                    total = 0
                    reads = {}

                # If the line has non-zero reads at the position,
                # write to the wig file, and update the total
                if count != 0:
                    reads[pos] = count
                    total += count



            # Write the final chromosome to the file
            print( current_chrom,total )
            out.write( "\n" + "variableStep chrom="+chrom )
            for my_pos,my_count in sorted(reads.items(),
                                        key=lambda x: x[0]):
                out.write("\n"+str(my_pos)+"\t"+str(my_count))



                # If the current line has a different chromosome,
                # write out a new header and update the current_chrom variable
                if chrom != current_chrom:
                    out.write( "\n" + "variableStep chrom="+chrom )
                    print( current_chrom,total )
                    current_chrom = chrom
                    total = 0


if __name__=="__main__":
    main()
