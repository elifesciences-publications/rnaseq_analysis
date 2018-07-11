import datetime as dt
import os
import pytest
import sys
import workflow
from modules.helpers import to_str

sys.path.append('.')


def test_process_config():
    output = workflow.process_config(config_file="tests/test_data/config")
    desired_output = {'bin_path': {'binbase': '/home/annasint/'}}
    assert output == desired_output


def test_submit_local_job():
    script = "java -jar /Users/annasintsova/tools/Trimmomatic-0.36/trimmomatic-0.36.jar " \
             "SE /Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq "\
             "/Users/annasintsova/git_repos/code/tests/test_data/UTI24_control_trimmed.fastq " \
             "ILLUMINACLIP:/Users/annasintsova/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true" \
             " SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0\n"

    workflow.submit_local_job(script)
    filename = "/Users/annasintsova/git_repos/code/tests/test_data/UTI24_control_trimmed.fastq"
    assert os.path.isfile(filename)
    assert os.path.getsize(filename) == 485388774


def test_run_trim_job(local_fastq_ref):
    # GIVEN:
    filename, _, today, config_dict, local = local_fastq_ref
    # WHEN:
    output_file_name = workflow.run_trim_job(filename, today, config_dict, local)[0]
    # THEN:
    assert os.path.isfile(output_file_name)
    assert os.path.getsize(output_file_name) != 0


def test_run_fastqc_job(local_fastq_ref, tmpdir):
    # GIVEN:
    filename, _, today, config_dict, local = local_fastq_ref
    # WHEN:
    actual_out_dir = workflow.run_fastqc_job(filename, today, config_dict, local)[0]
    # THEN:
    expected_filename = os.path.join(actual_out_dir,
                                     to_str(os.path.basename(filename).split(".fastq")[0]) + "_fastqc.html")
    assert os.path.join(str(tmpdir), "FastQC_results") == actual_out_dir
    assert len(os.listdir(actual_out_dir)) != 0
    assert os.path.isfile(expected_filename)


def test_run_build_index_job(local_fastq_ref):
    _, reference_genome, today, config_dict, local = local_fastq_ref
    bt2, _ = workflow.run_build_index_job(reference_genome,
                                          today, config_dict, local)
    assert os.path.isfile(bt2 + ".1.bt2")
    assert os.path.isfile(bt2 + ".4.bt2")


def test_run_align_job_local(local_fastq_ref):
    fastq_file, reference_genome, today, config_dict, local = local_fastq_ref
    bt2, _ = workflow.run_build_index_job(reference_genome, today, config_dict, local)
    sam_file, _ = workflow.run_alignment_job(fastq_file, bt2, config_dict, today, local)
    assert os.path.isfile(sam_file)
    assert os.path.getsize(sam_file) != 0


def test_run_sam_to_bam_conversion_and_sorting_local(local_sam):
    sam_file, today, config_dict, local = local_sam
    bam_file, _ = workflow.run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local)
    assert os.path.isfile(bam_file)
    assert os.path.isfile(bam_file + ".bai")


def test_run_count_job_bedtools_local(local_bam):
    bam, gff, today, config_dict, local = local_bam
    count_file = workflow.run_count_job_bedtools(gff, bam, config_dict, today, local)[0]
    assert os.path.isfile(count_file)
    assert os.path.getsize(count_file) != 0


###>>>>>>>

def test_run_counts_for_single_genome_bedtools_local(day): # for testing strand is False

    today = day
    bam_folder = "/Users/annasintsova/git_repos/code/data/alignments"
    gff = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
    config_dict = workflow.process_config("local_config")
    local = True
    workflow.run_counts_for_single_genome(gff, bam_folder, config_dict, today, local)
    file_names = workflow.find_files_in_a_tree(bam_folder, "bam")
    for file in file_names:
        count_file = file.split(".bam")[0] + "_counts_not_st.csv"
        assert os.path.getsize(count_file) != 0







def test_run_alignments_for_single_genome(local_mg1655_fastq_folder):
    genome, fastq_folder, today, config_dict, local = local_mg1655_fastq_folder
    workflow.run_alignments_for_single_genome(genome, fastq_folder, config_dict, today, local)

    file_names = workflow.find_files_in_a_tree(fastq_folder)  # todo test this function

    for file in file_names:
        bam_file = file.split(".")[0] + ".bam"
        assert os.path.isfile(bam_file)




if __name__ == "__main__":
    print("Hello!")
