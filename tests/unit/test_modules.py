import datetime as dt
import os
import sys
sys.path.append('.')

from modules import process_fastq
from modules import generate_PBS_script as pbs
from modules import align
from modules import samtools
from modules import counts


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


def test_generate_pbs(tmpdir):

    script_name = tmpdir.join("test.pbs")
    script = "mkdir TEST1"
    mem = '16gb'
    walltime ='24:00:00'
    pbs.generate_PBS_script(script_name, script, mem, walltime)
    assert os.path.isfile(script_name)
    with open(script_name, "r") as fo:
        text = fo.read()
    assert script in text
    assert "#PBS -l nodes=1:ppn=4,pmem={},walltime={}".format(mem, walltime) in text


def test_modules_fastqc(tmpdir):
    param_dict = {"bin": "fastqc_bin"}
    test_cases = ["not_fastq.txt", "zipped_file.fastq.gz", "fastq_file.fastq"]

    # Case 1
    script = process_fastq.fastqc(test_cases[0], tmpdir, param_dict)
    assert not script

    # Case 2
    script = process_fastq.fastqc(test_cases[1], tmpdir, param_dict)
    expected = "{} -o {} {} --extract\n".format(param_dict['bin'], tmpdir, test_cases[1])
    assert script == expected

    # Case 3
    script = process_fastq.fastqc(test_cases[2], tmpdir, param_dict)
    expected = "{} -o {} {}\n".format(param_dict['bin'], tmpdir, test_cases[2])
    assert script == expected


def test_modules_multiqc(tmpdir):
    input_dir = "input_dir"
    report = tmpdir.join(dt.datetime.today().strftime("%Y-%m-%d") + "_multiqc_report")
    param_dict = {"bin":"multiqc_bin"}
    script = process_fastq.multiqc(input_dir, tmpdir, param_dict)
    expected = "{} {} --force --filename {}\n".format(param_dict['bin'], input_dir, report)
    assert script == expected


def test_modules_align_bowtie():

    param_dict = {'bin': 'bowtie_bin', "seq_type": 'SE'}
    reference = "MG1655.fna"
    bt2_base = "MG1655_index"
    bin = param_dict['bin']
    script = align.build_bowtie_index(reference, bt2_base, param_dict)
    expected_script = "{}-build {} {}\n".format(bin, reference, bt2_base)
    assert script == expected_script

    fastq_file = "in.fastq"
    sam_file = "out.sam"
    script = align.bowtie_align(fastq_file, sam_file, bt2_base, param_dict)
    expected_script = "{} -x {} -U {} -S {}\n".format(bin, bt2_base, fastq_file, sam_file)
    assert script == expected_script

    param_dict['seq_type'] = 'PE'
    fastq_file = "in1.fastq in2.fastq"
    script = align.bowtie_align(fastq_file, sam_file, bt2_base, param_dict)
    expected_script = "{} -x {} -1 {} -2 {} -S {}\n".format(bin, bt2_base, fastq_file.split()[0],
                                                            fastq_file.split()[1], sam_file)
    assert script == expected_script


def test_modules_samtools():
    param_dict = {"bin":"samtools_bin"}
    sam_file = "in.sam"
    bam_file = "out.bam"
    sorted_name = "out_sorted.bam"
    script = samtools.sam2bam(sam_file, bam_file, param_dict)
    expected_script = "{0} view -b -o {1} -@ 4 {2}\n" \
                      "{0} sort -o {3} -@ 4 {1}\n" \
                      "{0} index -@ 4 {3}\n".format(param_dict['bin'], bam_file, sam_file, sorted_name)
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


def test_modules_counts_htseq_count():
    param_dict = {'bin':'htseq-count', 'form': 'bam', 'order': 'pos', 'attr': 'ID', 'mode': 'union',
                  'stranded': 'no', 'feature': 'CDS'}
    bin = param_dict['bin']
    form = param_dict["form"]
    order = param_dict["order"]
    attr = param_dict["attr"]
    mode = param_dict["mode"]
    stranded = param_dict["stranded"]
    feature = param_dict["feature"]
    bam = "input.bam"
    gff = 'reference'
    count_file = bam.split(".bam")[0] + "_counts.txt"
    script = counts.count_reads(bam, gff, param_dict)
    expected = "{0} -f {1} -r {2} -m {3}" \
               " -s {4} -t {5} -i {9} {6} {7} > {8}\n".format(bin, form, order, mode, stranded, feature,
                                                              bam, gff, count_file, attr)
    assert script == expected


def test_modules_counts_bedtool_count():
    param_dict = {'bin': "bedtools_bin", "strand": '-s'}
    bam_file = "input.bam"
    gff = "reference"
    count_file = bam_file.split(".bam")[0] + "_counts.txt"
    bedtools_bin = param_dict["bin"]
    strand = param_dict["strand"]
    expected = "{} coverage -a {} -b {} {}>{}".format(bedtools_bin, gff, bam_file, strand, count_file)
    script = counts.count_with_bedtools(gff, bam_file, param_dict)
    assert expected == script


# def test_modules_counts_count_with_bedtools_local():
#     gff = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
#     bam = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490_sorted.bam"
#     counts.count_with_bedtools_local(gff, bam, False, "locus_tag")
#
#     scrip(gff, bam, count_file, param_dict)
#     counts.count_with_bedtools_local(gff, bam, True, "locus_tag")
#     assert os.path.getsize(bam.split('.bam')[0] + "_counts_st.csv") != 0
#     assert os.path.getsize(bam.split('.bam')[0] + "_counts_not_st.csv") != 0