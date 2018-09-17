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


def count_reads(bam, gff, param_dict):
    """

    """
    bin = param_dict['bin']
    form = param_dict["form"]
    order = param_dict["order"]
    attr = param_dict["attr"]
    mode = param_dict["mode"]
    stranded = param_dict["stranded"]
    feature = param_dict["feature"]
    count_file = bam.split(".bam")[0] + "_counts.txt"
    script = "{0} -f {1} -r {2} -m {3}" \
             " -s {4} -t {5} -i {9} {6} {7} > {8}\n".format(bin, form, order, mode, stranded, feature,
                                                            bam, gff, count_file, attr)
    return script


def count_with_bedtools(gff, bam, param_dict):
    count_file = bam.split(".bam")[0] + "_counts.txt"
    bedtools_bin = param_dict["bin"]
    strand = param_dict["strand"]
    script = "{} coverage -a {} -b {} {}>{}".format(bedtools_bin, gff, bam, strand, count_file)
    return script


# def count_with_bedtools_local(gff, bam, count_file, param_dict):
#
#     strand = True if param_dict["strand"] == "-s" else False
#     feat = param_dict["feat"]
#     a = pybedtools.BedTool(gff)
#     b = pybedtools.BedTool(bam)
#     counts = a.coverage(b, counts=True, s=strand)
#     with open(count_file, "w") as fo:
#         for f in counts.features():
#             if feat not in str(f):
#                 continue
#             else:
#                 identifier = str(f[-2].split("{}=".format(feat))[1].split(";")[0].strip())
#                 fo.write("{},{}\n".format(identifier, f[-1]))
#
#     return count_file

if __name__ == "__main__":
    args = get_args().parse_args()
    param_dict = helpers.process_config(args.config)
    helpers.process_bedtools_count_output(args.input, args.output, param_dict)

