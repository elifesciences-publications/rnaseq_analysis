import pysam


def sam2bam(sam_file, bam_file,  param_dict):  # removed flagstat option
    samtools_bin = param_dict["bin"]
    sorted_name = bam_file.split(".")[0] + "_sorted.bam"
    script = "{0} view -b -o {1} -@ 4 {2}\n" \
             "{0} sort -o {3} -@ 4 {1}\n" \
             "{0} index -@ 4 {3}\n".format(samtools_bin, bam_file, sam_file, sorted_name)
    return script




# def flagstat(sam_file, param_dict):
#     script =
# todo finish flagstat option for flux
# Run flagstat, make summary, delete original flagstat files