import os


def build_bowtie_index(reference, bt2_base, param_dict):
    bowtie_bin = param_dict["bin"]
    script = "{}-build {} {}\n".format(bowtie_bin, reference, bt2_base)
    return script


def bowtie_align(fastq_file, sam_file, bt2_base, param_dict):
    assert 'fastq' in os.path.basename(fastq_file)
    bowtie_bin = param_dict["bin"]
    seq_type = param_dict["seq_type"]

    if seq_type == "PE":
        first = fastq_file.split()[0]
        second = fastq_file.split()[1]
        fastq_file = "-1 {} -2 {}".format(first, second)

    else:
        fastq_file = "-U {}".format(fastq_file)

    script = "{} -x {} {} -S {}\n".format(bowtie_bin, bt2_base, fastq_file, sam_file)
    return script
