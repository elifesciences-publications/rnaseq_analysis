
def sam2bam(sam_file, bam_file, flag_stat_file, samtools_bin):

    script = "{0} view -b -o {1} -@ 4 {2}\n" \
                "{0} flagstat {1} > {3}\n".format(samtools_bin, bam_file,
                                                  sam_file, flag_stat_file)

    sorted_name = bam_file.split(".")[0] + "_sorted.bam"
    script += "{0} sort -o {1} -@ 4 {2}\n" \
                      "{0} index -@ 4 {1}\n".format(samtools_bin, sorted_name, bam_file)
    return script



