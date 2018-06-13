import datetime as dt
import os

import subprocess
import sys
sys.path.append('.')
import workflow
from modules import process_fastq
from modules import generate_PBS_script as pbs
from modules import align
from modules.helpers import to_str


def test_modules_trimmomatic():
    input_file = "/Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq"
    output_directory = "/Users/annasintsova/git_repos/code/tests/test_data"
    suffix = os.path.basename(input_file).split(".fastq")[0]
    fastq_file_output = os.path.join(output_directory, suffix + "_trimmed.fastq")
    trimmomatic_bin = "/Users/annasintsova/tools/Trimmomatic-0.36/trimmomatic-0.36.jar"
    trimmomatic_adapters = "/Users/annasintsova/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa"

    script = process_fastq.trimmomatic(input_file, fastq_file_output,
                                       trimmomatic_bin, trimmomatic_adapters)
    assert script == "java -jar /Users/annasintsova/tools/Trimmomatic-0.36/trimmomatic-0.36.jar " \
                     "SE /Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq " \
                     "/Users/annasintsova/git_repos/code/tests/test_data/UTI24_control_trimmed.fastq " \
                     "ILLUMINACLIP:/Users/annasintsova/tools/" \
                     "Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true" \
                     " SLIDINGWINDOW:4:15 MINLEN:20 HEADCROP:0\n"



def test_generate_pbs():

    script_name = "/Users/annasintsova/git_repos/code/tests/test_data/test_script.pbs"
    script = "mkdir TEST1"

    pbs.generate_PBS_script(script_name, script)

    with open(script_name, "r") as fo:
        text = fo.read()
        print(text)
        assert script in text
        assert "#PBS -l nodes=1:ppn=4,pmem=4gb,walltime=12:00:00" in text


def test_modules_fastqc():


    # Case 1
    fastq_input_file = "/Users/annasintsova/git_repos/code/data/reads/UTI24_control.txt"
    out_dir = "/Users/annasintsova/git_repos/code/tests/test_data"
    fastqc_bin = "/Users/annasintsova/tools/FastQC/fastqc"
    script = process_fastq.fastqc(fastq_input_file, out_dir, fastqc_bin)
    assert not script
    # Case 2
    fastq_input_file = "/Users/annasintsova/git_repos/code/data/reads/UTI24_control_copy.fastq.gz"
    script = process_fastq.fastqc(fastq_input_file, out_dir, fastqc_bin)
    expected = "/Users/annasintsova/tools/FastQC/fastqc -o " \
               "/Users/annasintsova/git_repos/code/tests/test_data " \
               "/Users/annasintsova/git_repos/code/data/reads/UTI24_control_copy.fastq.gz --extract\n"
    assert script == expected
    # Case 3
    fastq_input_file = "/Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq"
    script = process_fastq.fastqc(fastq_input_file, out_dir, fastqc_bin)
    expected = "/Users/annasintsova/tools/FastQC/fastqc -o " \
               "/Users/annasintsova/git_repos/code/tests/test_data " \
               "/Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq\n"
    assert script == expected

def test_modules_align_build_index():
    reference = "data/ref/MG1655.fna"
    bt2_base = "data/ref/MG1655_index"
    bowtie_bin = ""
    script = align.build_bowtie_index(reference, bt2_base,
                                      bowtie_bin)
    expected_script = "bowtie2-build data/ref/MG1655.fna data/ref/MG1655_index\n"
    assert script == expected_script

def test_modules_align():
    fastq_file = "/Users/annasintsova/git_repos/code/data/reads/HM86_UR_copy.fastq.gz"
    file_name = "/Users/annasintsova/git_repos/code/data/reads/HM86_UR_copy.fastq"
    sam_file_name = "/Users/annasintsova/git_repos/code/data/reads/HM86_UR_copy.sam"
    bt2_base = "/Users/annasintsova/git_repos/code/data/ref/MG1655_index"
    bowtie_bin = ""
    script = align.bowtie_align(fastq_file,
                                sam_file_name,
                                bt2_base, bowtie_bin)
    expected_script = "gunzip {}\n{}bowtie2 -x {} "\
                      "-U {} -S {}\n".format(fastq_file, bowtie_bin, bt2_base,
                                             file_name, sam_file_name)
    assert script == expected_script
