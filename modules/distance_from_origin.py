# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:11:11 2017

@author: annasintsova
"""

# Only working with genes etc locationd on first contig, eventually hoping to have a close chromosome
import pandas as pd
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


genome_info = {'HM1':(5284090,3631354),#length, origin
'HM3':(4729556,2181548),
'HM6':(5344257,2822711),
'HM7':(4886521,3306317),
'HM14':(5037986,3831641),
'HM17':(5087802,2557959),
'HM43':(4901204,2345263),
'HM54':(5546219,788250),
'HM56':(5040698,761996),
'HM57':(5152486,2827113),
'HM66':(5168961,2648967),
'HM68':(5054821,3298219),
'HM86':(5220431,686494)}




# differential_expression_file
# crossRef file
# genome annotation




#given PROKKA calculate distance to origin

#for given annotation file, read in genes that are on first contig and save their start point
# for each calculate distance to the origin

def calculate_distance_for_genome(genome, gff_file, genome_info = genome_info):
    length = genome_info[genome][0]
    origin = genome_info[genome][1]
    gene_distance = {}
    with open(gff_file, "r") as fh:
        for line in fh:
            if line.startswith("{}_1".format(genome)):
                gene_id = line.split("\t")[8].split(";")[0].strip("ID=")
                gene_position = int(line.split("\t")[3])
                distance = distance_to_origin(gene_position, origin, length)
                gene_distance[gene_id] = distance
    return gene_distance


def get_differentially_expressed_genes(deseq_file, FC = 0):
    de_genes = {}
    with open(deseq_file, "r") as fh:
        fh.readline()
        for line in fh:
            gene = line.split(",")[0]
            lfc = float(line.split(",")[2])
            if not FC:
                de_genes[gene] = lfc
            else:
                if lfc > FC or lfc < FC*(-1):
                    de_genes[gene] = lfc
    return de_genes


def get_prokkas(genome, cross_ref_file):

    df = pd.read_csv(cross_ref_file, index_col=0)
    df = df.dropna(subset=["MG1655"])
    df.set_index("MG1655", inplace=True)
    genes = df[genome]
    bnum_to_prokkas = genes.to_dict()
    return bnum_to_prokkas


def distance_to_origin(gene_position, origin, length):
    if gene_position < origin:
        distance = 1 + ((gene_position - origin)/length)
        #distance2 = (length - (origin - gene_position))/length
        return round(distance, 2)
    else:
        distance = (gene_position - origin)/length
        return round(distance, 2)

def FC_vs_distance(de_genes, gene_distance, get_prokkas, out_file = ''):


    distances = []
    names = []
    LFCs = []
    for key, val in de_genes.items():
        prokka = get_prokkas[key]
        if prokka in gene_distance.keys():
            distance = gene_distance[prokka]
            names.append(prokka)
            distances.append(distance)
            LFCs.append(val)
    df = pd.DataFrame({"Distance":distances,
                       "LFC": LFCs},
                       index = names)
    df.to_csv(out_file)
    print(df[0:3])

def calculate_fc_vs_distance(genome, gff_file, cross_ref_file, out_file = ''):

    prokkas = get_prokkas(genome, cross_ref_file)
    de_genes = get_differentially_expressed_genes(deseq_file, FC = 1)

    gene_distance = calculate_distance_for_genome(genome, gff_file, genome_info = genome_info)
    print(gene_distance)
    FC_vs_distance(de_genes, gene_distance, prokkas, out_file)

if __name__ == "__main__":
    print(distance_to_origin(163441,  3306317, 4886521))

    deseq_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/analysis/DE/2018-03-06-DEseq2-analysis/data/best_strains_DEseq.csv"

    hm86_gff = "/Users/annasintsova/Downloads/RNA_reference_genomes/HM56.gff"

    cross_ref_file = "/Users/annasintsova/git_repos/HUTI-RNAseq/data/get_homologs_output/C50_S90_e0_/run_C50_S90_e0__pan_C50_S90/2018-02-26_pangenome_matrix_t0_crossRef.csv"

    calculate_fc_vs_distance("HM56", hm86_gff, cross_ref_file, "tests/hm56_distances.csv")



