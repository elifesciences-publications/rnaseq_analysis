import os


def build_bowtie_index(reference, bt2_base, bowtie_bin):
    script = "{}bowtie2-build {} {}\n".format(bowtie_bin, reference, bt2_base)
    return script


def bowtie_align(fastq_file,
                 sam_file_name,
                 bt2_base, bowtie_bin):
    script = ""
    assert 'fastq' in os.path.basename(fastq_file)
    if fastq_file.endswith(".gz"):
        script += "gunzip {}\n".format(fastq_file)
        file_name = fastq_file.split(".gz")[0]
    else:
        file_name = fastq_file
    script += "{}bowtie2 -x {} {} -S {}\n".format(bowtie_bin, bt2_base,
                                                  file_name, sam_file_name)
    return script
