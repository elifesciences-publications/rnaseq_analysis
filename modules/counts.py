import pybedtools
import os


def count_reads(bam, gff, config_dict):
    """

    """
    htseq_path = config_dict["HTSeq"]["bin"]
    form = config_dict["HTSeq"]["form"]
    order = config_dict["HTSeq"]["order"]
    attr = config_dict["HTSeq"]["attr"]
    mode = config_dict["HTSeq"]["mode"]
    stranded = config_dict["HTSeq"]["stranded"]
    feature = config_dict["HTSeq"]["feature"] # todo update flux config
    count_file = bam.split(".bam")[0] + "_counts.txt"
    script = "htseq-count -f {0} -r {1} -m {2}" \
             " -s {3} -t {4} -i {8} {5} {6} > {7}\n".format(form, order,
                                                            mode, stranded, feature, bam,
                                                            gff, count_file, attr)
    return script


def count_with_bedtools_flux(gff, bam, count_file, config_dict, strand=False):

    bedtools_path = config_dict["bedtools"]["bin"]
    s = " -s" if strand else ""
    script = "{} coverage -a {} -b {} -counts{}>{}".format(bedtools_path, gff, bam, s, count_file)
    return script


def count_with_bedtools_local(gff, bam, count_file, strand=False, feat="locus_tag"):
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

# if __name__ == "__main__":
#     gff = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
#     bam = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490_sorted.bam"
#     count_with_bedtools_local(gff, bam, strand=False, feat="locus_tag")

