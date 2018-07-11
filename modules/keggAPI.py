"""
__author__ = annasint

for each gene see if there's NAME, PATHWAY, BRITE field

ENTRY
NAME
DEFINITION
ORTHOLOGY
ORGANISM
BRITE
POSITION
MOTIF
DBLINKS
STRUCTURE
AASEQ
NTSEQ
///


"""
import urllib.request
from collections import OrderedDict
import argparse
import datetime as dt
import os


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--gene_list_file',  help='file with list of genes',
                        type=str, required=True)
    parser.add_argument('-g', '--genome', help='Genome', type=str, required=True)
    parser.add_argument('-o', '--output_prefix', type=str, required=True)
    return parser


def get_genes(genome, gene_list_txt, sep=",", hdr=True, id_col=0):
    """
    :param genome: string, KEGG genome designation, ex. ecc for CFT073, eco for MG1655
    :param gene_list_txt: list of genes for which info/seq is needed, have to Gene entry name or accession (from KEGG API)
    :param sep: how fields are separated in the gene_list_txt
    :param hdr: whether gene_list_txt has a header
    :param id_col: where gene name/tag located in the file
    :return: list KEGG API compatible gene requests, ex. genome:geneid
    """

    with open(gene_list_txt) as gl:
        if hdr:
            gl.readline() # skipping the header
        try:
            genes = [line.rstrip().split(sep)[id_col] for line in gl]
        except:
            return []
    return ["{}:{}".format(genome, gene) for gene in genes]


def getgeneInfo(gene_call_list, output_prefix):

    all_gene_info = []
    if not len(gene_call_list):
        return "No genes in the list"
    for gene in gene_call_list:
        try:
            url = "http://rest.kegg.jp/get/{}".format(gene)
            with urllib.request.urlopen(url) as f:
                lines = f.read().decode('utf-8').splitlines()

            gene_info = OrderedDict([('locus_tag', ''),
                                     ('name', ''),
                                     ('definition', ''),
                                     ('pathway', ''),
                                     ('aa_seq', ''),
                                     ('nt_seq', '')])  # can add more if needed

            for i in range(len(lines)):
                field = lines[i].split()[0]
                if field == 'ENTRY':
                    gene_info['locus_tag'] = lines[i].split()[1]

                elif field == "NAME":
                    p = lines[i].split("NAME        ")[1]
                    gene_info['name'] = p

                elif field == "DEFINITION":

                    p = lines[i].split("DEFINITION  ")[1]
                    print(p)
                    gene_info['definition'] = p

                elif field == "PATHWAY":
                    p = lines[i].split("PATHWAY     ")[1]
                    count = 1
                    while lines[i+count].startswith("           ") and len(lines[i+count]) > 0:
                        p += lines[i+count]
                        count += 1
                    gene_info['pathway'] = ", ".join(p.split("            ")).rstrip(', ')

                elif field == "AASEQ":
                    p = ''
                    count = 1
                    while lines[i+count].startswith("           ") and len(lines[i+count]) > 0:
                        p += lines[i+count]
                        count +=1
                    gene_info['aa_seq'] = p.replace("            ", '')

                elif field == "NTSEQ":
                    p = ''
                    count = 1
                    while lines[i+count].startswith("           ") and len(lines[i+count]) > 0:
                        p += lines[i+count]
                        count +=1
                    gene_info['nt_seq'] = p.replace("            ", '')


            info_str = "\t".join([v for k,v in gene_info.items() if k in ["locus_tag", "name", "definition","pathway"]]).rstrip('\t')
            all_gene_info.append(gene_info)
            aa_str = ">{}\n{}".format(gene, gene_info['aa_seq'])
            nt_str = ">{}\n{}".format(gene, gene_info['nt_seq'])
            # Want to write to file right away, in case connection breaks
            with open(output_prefix+"_info.tab", "a") as of: # todo refactor this
                of.write(info_str+"\n")
            with open(output_prefix+"_aa.fasta", "a") as of:
                of.write(aa_str+"\n")
            with open(output_prefix+"_nt.fasta", "a") as of:
                of.write(nt_str+"\n")
        except:
            with open(output_prefix+"_info.tab", "a") as of:
                of.write("{}\tNaN\tNaN\tNaN\n".format(gene.split(":")[1]))

            continue
    return all_gene_info

def searchKEGGGenome(genome="eco", search_term="ribosomal subunit", out_dir="."):
    url = "http://rest.kegg.jp/list/{}".format(genome)
    search_results = []
    with urllib.request.urlopen(url) as f:
        lines = f.read().decode('utf-8').splitlines()
    for line in lines:
        if search_term in line:
            search_results.append(line)
    out_file = "{}_{}_{}.kegg.results".format(dt.datetime.today().strftime("%Y-%m-%d"),
                                              genome, search_term.replace(" ", "_"))
    with open(os.path.join(out_dir, out_file), "w+") as fo:
        for r in search_results:
            fo.write(r + "\n")
    return [line.split()[0].split(":")[1] for line in search_results]


def explore_kegg():
    url = "http://rest.kegg.jp/list/ecc"
    with urllib.request.urlopen(url) as f:
        lines = f.read().decode('utf-8').splitlines()
        genes = [line.split()[0] for line in lines]

    print("Number of genes: {}".format(len(genes)))
    pathway_field = 0
    brite_field = 0
    name_field = 0
    definition_field = 0
    organism_field = 0
    aa_seq_field = 0

    info = []
    for i in range(0, len(genes)+1, 10):
        gene_request = "+".join(genes[i:i+10])
        url = "http://rest.kegg.jp/get/{}".format(gene_request)
        with urllib.request.urlopen(url) as f:
            lines = f.read().decode('utf-8').splitlines()
            info.append(lines)
    for lines in info:

        for l in lines:
            print(l)
            if l.split()[0] == "NAME":
                name_field += 1
            elif l.split()[0] == "DEFINITION":
                definition_field += 1
            elif l.split()[0] == "ORGANISM":
                organism_field += 1
            elif l.split()[0] == "PATHWAY":
                pathway_field += 1
            elif l.split()[0] == "BRITE":
                brite_field += 1
            elif l.split()[0] == "AASEQ":
                aa_seq_field += 1
    print("Number of NAME fields: {}".format(name_field))
    print("Number of DEFINITION fields: {}".format(definition_field))
    print("Number of ORGANISM fields: {}".format(organism_field))
    print("Number of PATHWAY fields: {}".format(pathway_field))
    print("Number of BRITE fields: {}".format(brite_field))
    print("Number of AASEQ fields: {}".format(aa_seq_field))


def process_gene_id_file(genome, gene_list_txt, output_prefix, sep=",", hdr=True, id_col=0):
    gene_list = get_genes(genome, gene_list_txt, sep, hdr, id_col)
    info = getgeneInfo(gene_list, output_prefix)
    return info


if __name__ == "__main__":

    args = parser().parse_args()
    gene_call_list = get_genes(args.genome, args.gene_list_file)
    getgeneInfo(gene_call_list, args.output_prefix)
