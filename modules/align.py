import os


def build_bowtie_index(reference, bt2_base, output_directory,
                       bowtie_bin):
    script = "cd {}\n".format(output_directory)
    script += "{}/bowtie2-build {} {}\n".format(bowtie_bin, reference, bt2_base)
    return script


def bowtie_align(fastq_file,  output_directory,
                 sam_file_name,
                 bt2_base, bowtie_bin):
    script = "cd {}\n".format(output_directory)
    if 'fastq' not in os.path.basename(fastq_file):
        return "Skipping {}, not a fastq file".format(fastq_file)
    elif fastq_file.endswith(".gz"):
        script += "gunzip {}\n".format(fastq_file)
        file_name = fastq_file.split(".gz")[0]
    else:
        file_name = fastq_file
    script += "{}/bowtie2 -x {} -U {} -S {}\n".format(bowtie_bin, bt2_base,
                                                      file_name, sam_file_name)
    return script
