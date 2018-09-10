import datetime as dt
import os

import subprocess
import sys
sys.path.append('.')
import workflow
from modules import process_fastq
from modules import generate_PBS_script as pbs
from modules import align
from modules import samtools
from modules import counts
from modules.helpers import to_str


def test_modules_trimmomatic():
    # SE
    input_file = "input.fastq"
    output_file = "output.fastq"
    param_dict = {"seq_type": "SE",
                  "bin": "trimmomatic-0.36.jar",
                  "adapters_se": "TruSeq3-SE.fa",
                  "adapters_pe": "TruSeq3-PE-2.fa",
                  "headcrop": "20",
                  "crop": "0",
                  "seed_mismatches": "2",
                  "palindrome_clipthreshold": "30",
                  "simple_clipthreshold": "10",
                  "minadapterlength": "8",
                  "keep_both_reads": "true",
                  "window_size": "4",
                  "window_size_quality": "20",
                  "minlength": "40"}

    script = process_fastq.trimmomatic(input_file, output_file, param_dict)
    expected_script = "java -jar {bin} {seq_type} {0} {1} ILLUMINACLIP:{adapters_se}:" \
                      "{seed_mismatches}:{palindrome_clipthreshold}:{simple_clipthreshold} " \
                      "SLIDINGWINDOW:{window_size}:{window_size_quality} MINLEN:{minlength} " \
                      "HEADCROP:{headcrop}\n".format(input_file, output_file, **param_dict)

    assert script == expected_script
    # PE
    param_dict["seq_type"] = "PE"
    script = process_fastq.trimmomatic(input_file, output_file, param_dict)
    expected_script = "java -jar {bin} {seq_type} {0} {1} ILLUMINACLIP:{adapters_pe}:" \
                      "{seed_mismatches}:{palindrome_clipthreshold}:{simple_clipthreshold}:" \
                      "{minadapterlength}:{keep_both_reads} " \
                      "SLIDINGWINDOW:{window_size}:{window_size_quality} MINLEN:{minlength} " \
                      "HEADCROP:{headcrop}\n".format(input_file, output_file, **param_dict)

    assert script == expected_script


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

def test_modules_samtools_get_bam_stats():
    bam = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490_sorted.bam"
    total = 1605460
    mapped = 1347672
    pcnt_mapped = 83.94
    actual_total, actual_mapped, actual_pcnt_mapped = samtools.get_bam_stats(bam)
    assert total == actual_total
    assert mapped == actual_mapped
    assert pcnt_mapped == actual_pcnt_mapped




def test_modules_counts_count_with_bedtools_local():
    gff = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
    bam = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490_sorted.bam"
    counts.count_with_bedtools_local(gff, bam, False, "locus_tag")
    counts.count_with_bedtools_local(gff, bam, True, "locus_tag")
    assert os.path.getsize(bam.split('.bam')[0] + "_counts_st.csv") != 0
    assert os.path.getsize(bam.split('.bam')[0] + "_counts_not_st.csv") != 0