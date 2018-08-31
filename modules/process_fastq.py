import datetime as dt
import os

# FASTQC, TRIMMOMATIC

# All of these need to be adjusted to accomodate PE situation


def fastqc(fastq_input_file,  out_dir, fastqc_bin):
    """

    """
    if "fastq" not in os.path.basename(fastq_input_file):
        print("Not a fastq file")
        return ''
    elif ".gz" in fastq_input_file:
        return "{} -o {} {} --extract\n".format(fastqc_bin, out_dir, fastq_input_file)
    else:
        return "{} -o {} {}\n".format(fastqc_bin, out_dir, fastq_input_file)


def multiqc(input_directory, output_directory=''):

    if not output_directory:
        output_directory = input_directory
    report_name = os.path.join(output_directory,
                               (dt.datetime.now().strftime("%Y-%m-%d") + "_multiqc_report"))

    return "multiqc {} --force --filename {}\n".format(input_directory, report_name)
    # todo multiqc is not installed in the bin


def trimmomatic(fastq_file_input, fastq_file_output, config_dict):

        # if PE fastq_file_input a string of forward and reverse file names
        # seperated by a whitespace
        # same for the output

        trimmomatic_bin = config_dict["Trimmomatic"]["bin"]
        trimmomatic_adapters = config_dict["Trimmomatic"]["adapters"]
        seq_type = config_dict["sequencing"]["type"]
        headcrop = config_dict["Trimmomatic"]["headcrop"]
        minlen = config_dict["Trimmomatic"]["minlen"]
        return "java -jar {} " \
               "{} {} {} ILLUMINACLIP:{}:2:30:10:8:true" \
               " SLIDINGWINDOW:4:15 MINLEN:{} HEADCROP:{}\n".format(trimmomatic_bin,
                                                                    seq_type,
                                                                    fastq_file_input,
                                                                    fastq_file_output,
                                                                    trimmomatic_adapters,
                                                                    minlen,
                                                                    headcrop)

