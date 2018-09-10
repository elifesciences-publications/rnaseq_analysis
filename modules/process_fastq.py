"""

Functions to run FastQC, MultiQC, and Trimmomatic

"""
import datetime as dt
import os


def fastqc(fastq_input_file, out_dir, param_dict):
    """
    Takes in fastq_file, output directory and parameter dictionary
    Returns command to run FastQC on that fastq file

    """
    if "fastq" not in os.path.basename(fastq_input_file):
        print("Not a fastq file")
        return ''
    elif ".gz" in fastq_input_file:
        return "{} -o {} {} --extract\n".format(param_dict["bin"], out_dir, fastq_input_file)
    else:
        return "{} -o {} {}\n".format(param_dict["bin"], out_dir, fastq_input_file)


def multiqc(input_directory, output_directory, param_dict):

    report_name = os.path.join(output_directory,
                               (dt.datetime.now().strftime("%Y-%m-%d") + "_multiqc_report"))

    return "{} {} --force --filename {}\n".format(param_dict["bin"], input_directory, report_name)
    # todo multiqc is not installed in the bin


def trimmomatic(fastq_file_input, fastq_file_output, param_dict):  # todo refactor further
    """Takes in input, output, dictionary of optional parameters
       Returns script to run trimmomatic
       If PE: fastq_file_input a string of forward and reverse file names seperated by a whitespace
       Same for the output

    """
    seq_type = param_dict["seq_type"]
    path_string = 'java -jar {bin} {seq_type}'.format(**param_dict)
    sliding_string = 'SLIDINGWINDOW:{window_size}:{window_size_quality}'.format(**param_dict)
    minlen_string = 'MINLEN:{minlength}'.format(**param_dict)

    if seq_type == "SE":
        illumina_string = 'ILLUMINACLIP:{adapters_se}:{seed_mismatches}:{palindrome_clipthreshold}:' \
                          '{simple_clipthreshold}'.format(**param_dict)
    else:
        illumina_string = 'ILLUMINACLIP:{adapters_pe}:{seed_mismatches}:{palindrome_clipthreshold}:' \
                          '{simple_clipthreshold}:{minadapterlength}:{keep_both_reads}'.format(**param_dict)

    headcrop = param_dict["headcrop"]
    crop = param_dict["crop"]
    if int(crop):
        headcrop_string = ""
        crop_string = 'CROP:{crop} '.format(**param_dict)
    elif int(headcrop):
        headcrop_string = ' HEADCROP:{headcrop}'.format(**param_dict)
        crop_string = ""
    else:
        headcrop_string = ""
        crop_string = ""

    return "{} {} {} {}{} {} {}{}\n".format(path_string, fastq_file_input, fastq_file_output,
                                            crop_string, illumina_string, sliding_string,
                                            minlen_string, headcrop_string)
