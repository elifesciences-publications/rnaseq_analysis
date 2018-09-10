import argparse
import pybedtools
import sys
sys.path.append("/Users/annasintsova/git_repos/code/modules/")
import helpers


def get_args():
    parser = argparse.ArgumentParser("RNASeq Pipeline\n")
    parser.add_argument("-i", "--input",  required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-c", "--config", required=True)
    return parser


def count_reads(bam, gff, config_dict):
    """

    """
    form = config_dict["HTSeq"]["form"]
    order = config_dict["HTSeq"]["order"]
    attr = config_dict["HTSeq"]["attr"]
    mode = config_dict["HTSeq"]["mode"]
    stranded = config_dict["HTSeq"]["stranded"]
    feature = config_dict["HTSeq"]["feature"]  # todo update flux config
    count_file = bam.split(".bam")[0] + "_counts.txt"
    script = "htseq-count -f {0} -r {1} -m {2}" \
             " -s {3} -t {4} -i {8} {5} {6} > {7}\n".format(form, order,
                                                            mode, stranded, feature, bam,
                                                            gff, count_file, attr)
    return script


def count_with_bedtools(gff, bam, count_file, param_dict):

    bedtools_bin = param_dict["bin"]
    strand = param_dict["strand"]
    script = "{} coverage -a {} -b {} {}>{}".format(bedtools_bin, gff, bam, strand, count_file)
    return script


def count_with_bedtools_local(gff, bam, count_file, param_dict):

    strand = True if param_dict["strand"] == "-s" else False
    feat = param_dict["feat"]
    a = pybedtools.BedTool(gff)
    b = pybedtools.BedTool(bam)
    counts = a.coverage(b, counts=True, s=strand)
    with open(count_file, "w") as fo:
        for f in counts.features():
            if feat not in str(f):
                continue
            else:
                identifier = str(f[-2].split("{}=".format(feat))[1].split(";")[0].strip())
                fo.write("{},{}\n".format(identifier, f[-1]))

    return count_file


def process_bedtools_count_output(count_file, count_file_out, param_dict):
    feat = param_dict["feat"]
    count_dict = {}
    with open(count_file, "r") as fh:
        for line in fh:
            if feat not in line:
                continue
            else:
                gene_info = line.split("\t")[8].split(feat+"=")[1].split(";")[0].strip()
                counts = line.split("\t")[9].strip()
                coverage = line.split("\t")[-1].strip()
                count_dict[gene_info] = (counts, coverage)
    with open(count_file_out, "w") as fo:
        for key, val in count_dict.items():
            fo.write("{},{},{}\n".format(key, val[0], val[1]))
    return count_file_out


if __name__ == "__main__":
    args = get_args().parse_args()
    param_dict = helpers.process_config(args.config)
    process_bedtools_count_output(args.input, args.output, param_dict)

