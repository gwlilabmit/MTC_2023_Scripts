import sys
import pandas as pd
import os

def get_dataframe(gen, cds_dict, genomes, file):

    for g in genomes:
        if g in gen:
            cds_df = cds_dict[g].copy()
    if 'minus' in file:
        cds_df = cds_df[cds_df.Strand == '+']
    else:
        cds_df = cds_df[cds_df.Strand == '-']
    cds_df = cds_df.reset_index()
    cds_df = cds_df.drop(columns = ['index'])
    return cds_df

def gene_count(wig_prefix, cds_dict, names, genomes, save_prefix):

    print('Aligning for sample {}'.format(wig_prefix))
    file_names = [wig_prefix + '_plus.wig', wig_prefix + '_minus.wig']
    dfs =[]
    for file in file_names:
        print('Now adding reads from file {}'.format(file))
        with open(file, 'r') as f:
            l = f.readline()
            l = f.readline()
            l = l.rstrip()
            gen = l[l.rfind('=')+1:]
            temp_df = get_dataframe(gen, cds_dict, genomes, file)
            num_genes = len(temp_df)
            cur_ind = 0
            for l in f:
                l = l.rstrip()
                if '=' in l:
                    dfs.append(temp_df)
                    gen = l[l.rfind('=')+1:]
                    temp_df = get_dataframe(gen, cds_dict, genomes, file)
                    print(temp_df.head)
                    num_genes = len(temp_df)
                    cur_ind = 0
                else:
                    pos, counts = l.split('\t')
                    pos = int(float((pos)))
                    counts = int(float(counts))
                    while cur_ind <  num_genes and pos > temp_df.iloc[cur_ind].Stop:
                        cur_ind += 1
                    if cur_ind < num_genes:
                        if pos > temp_df.iloc[cur_ind].Start:
                            temp_df.at[cur_ind, 'Counts'] += counts
        dfs.append(temp_df)
    merge_df = pd.concat(dfs, ignore_index=True, sort=False)
    merge_df.to_csv(''.join([save_prefix, '_dataframe.txt']))


def make_CDS_df(cds_file, name):
    """
    Makes a dataframe for each CDS file for easy access downstream.
        -cds_file = str with the path of the cds_file of interest.

    Returns a pandas dataframe with the relevant columsn of the original CDS.
    """

    cds_saved_file = ''.join([cds_file[:cds_file.rfind('.txt')], '_processed.txt'])
    if os.path.exists(cds_saved_file):
        df = pd.read_csv(cds_saved_file, index_col = 0)
        return df

    cds_df = pd.DataFrame(columns = ['Name', 'Strand', 'Start', 'Stop', 'Note'])

    with open(cds_file, 'r') as cds:
        cds.readline()
        cds.readline()

        for l in cds:
            space_break = l.split('\t')
            if len(space_break) > 3 and space_break[2] in ['CDS', 'rRNA']:
                start = int(space_break[3])
                end = int(space_break[4])
                strand = space_break[6]
                colon_break = space_break[8].split(';')
                gene_name = 'error'
                note = ''
                for i in colon_break:
                    if i.startswith('gene='):
                        gene_name = i[i.index('=')+1:]
                    if i.startswith('Note='):
                        note = i[i.index('=')+1:]
                    if i.startswith('product='):
                        product = i[i.index('=')+1:]
                if note == '':
                    note = product
                cds_df = cds_df.append({'Name' : gene_name, 'Strand' : strand, 'Start' : start, 'Stop': end, 'Note': note},
                        ignore_index = True)

        cds_df = cds_df.sort_values(by = ['Strand', 'Start'])
        cds_df = cds_df.reset_index()
        cds_df = cds_df.drop(columns = ['index'])
        cds_df['Genome'] = [name for i in range(len(cds_df.Name))]
        cds_df['Counts'] = [0 for _ in range(len(cds_df.Name))]
        cds_df.to_csv(cds_saved_file)

    return cds_df


if __name__ == '__main__':
    """
    Program starts here when called from command line.
    Expected arguments (in order):
        -wig_file = string of where name of current wig file
        -cds_folder = string path where CDS files are kept
        -cds_files = `,` seperated str with names of CDS files
        -genomes = `,` seperated str with the corresponding genomes of the CDS files (eg NC_000913.2)
        -names = `,` seperated str with names for each genome/CDS pair (eg ecoli)
        -prefix = str with prefix for save location of the dataframe file.
    """

    wig_prefix = sys.argv[1]
    cds_folder = sys.argv[2]
    cds_files = sys.argv[3].split(",")
    genomes = sys.argv[4].split(",")
    names = sys.argv[5].split(",")
    save_prefix = sys.argv[6]

    cds_dict = {
        key: make_CDS_df(''.join([cds_folder, cds_files[ind]]), genomes[ind])
        for ind, key in enumerate(genomes)
    }
    gene_count(wig_prefix, cds_dict, names, genomes, save_prefix)
