import pysam
#import pysamstats
#import matplotlib.pyplot as plt


def sam2bam(sam_file, bam_file,  samtools_bin):  # removed flagstat option

    script = "{0} view -b -o {1} -@ 4 {2}\n".format(samtools_bin, bam_file, sam_file)
    sorted_name = bam_file.split(".")[0] + "_sorted.bam"
    script += "{0} sort -o {1} -@ 4 {2}\n" \
                      "{0} index -@ 4 {1}\n".format(samtools_bin, sorted_name, bam_file)
    return script


def get_bam_stats(bam):
    mybam = pysam.AlignmentFile(bam)
    mapped = mybam.mapped
    total = mybam.mapped + mybam.unmapped
    pnt_mapped = mapped/total*100
    return total, mapped, round(pnt_mapped, 2)



# if __name__ == "__main__":
