import datetime as dt
import os
import pytest
import sys
import workflow
from modules import helpers

sys.path.append('.')


def test_get_second():
    first = "first_1.fastq"
    second = helpers.get_second(first)
    expected_second = "first_2.fastq"
    assert second == expected_second
    first = "first.forward.fastq"
    second = helpers.get_second(first)
    expected_second = "first.reverse.fastq"
    assert second == expected_second
    first = "first_2.fastq"
    second = helpers.get_second(first)
    expected_second = "first_2.fastq"
    assert second == expected_second
    first = "/Users/some/path/to/dir/first_1.fastq"
    second = helpers.get_second(first)
    expected_second = "/Users/some/path/to/dir/first_2.fastq"
    assert second == expected_second


def test_process_config():
    output = helpers.process_config(config_file="tests/test_config")
    actual_output = output['bin_path']['binbase']
    expected_output = '/home/annasint/'
    assert actual_output == expected_output


def test_submit_local_job():
    script = "echo Hello\necho World!\n"
    output = workflow.submit_local_job(script)
    expected_output = ["Hello", "World!"]
    assert output == expected_output


def test_run_trim_job(local_fastq_ref):
    # GIVEN:
    filename, _, config_dict, local = local_fastq_ref
    # Because of the test fastq file, this will only pass
    # if Trimmomatic MINLEN is set to 20
    # And no cropping, otherwise second assertion fails

    # WHEN:
    output_file_name = workflow.run_trim_job(filename, config_dict, local)[0]
    # THEN:
    assert os.path.isfile(output_file_name)
    assert os.path.getsize(output_file_name) != 0
# todo test run_trim_job for PE


def test_run_fastqc_job(local_fastq_ref, tmpdir):
    # GIVEN:
    filename, _, config_dict, local = local_fastq_ref
    sample_id = helpers.to_str(os.path.basename(filename).split(".fastq")[0])
    # WHEN:
    actual_out_dir = workflow.run_fastqc_job(filename, config_dict, local)[0]
    # THEN:
    expected_filename = os.path.join(actual_out_dir, sample_id + "_fastqc.html")
    assert os.path.isdir(actual_out_dir)
    assert len(os.listdir(actual_out_dir)) != 0
    assert os.path.isfile(expected_filename)

# todo test run_fastqc_job for PE

# todo test multiqc


def test_run_build_index_job(local_fastq_ref):
    _, ref_genome, config_dict, local = local_fastq_ref
    index_path, _ = workflow.run_build_index_job(ref_genome, config_dict, local)
    assert os.path.isfile(index_path + ".1.bt2")
    assert os.path.isfile(index_path + ".4.bt2")


def test_run_align_job(local_fastq_ref):
    fastq_file, ref_genome, config_dict, local = local_fastq_ref
    bt2, _ = workflow.run_build_index_job(ref_genome, config_dict, local)
    sam_file, _ = workflow.run_alignment_job(fastq_file, bt2, config_dict, local)
    assert os.path.isfile(sam_file)
    assert os.path.getsize(sam_file) != 0

# todo test run_align for PE samples


def test_run_sam_to_bam_conversion_and_sorting(local_sam):
    sam_file, today, config_dict, local = local_sam
    bam_file, _ = workflow.run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local)
    assert os.path.isfile(bam_file)
    assert os.path.isfile(bam_file + ".bai")


# def test_run_count_job_bedtools_local(local_bam):
#     bam, gff, today, config_dict, local = local_bam
#     count_file = workflow.run_count_job_bedtools(gff, bam, config_dict, today, local)[0]
#     assert os.path.isfile(count_file)
#     assert os.path.getsize(count_file) != 0


###>>>>>>>
#
# def test_run_counts_for_single_genome_bedtools_local(day): # for testing strand is False
#
#     today = day
#     bam_folder = "/Users/annasintsova/git_repos/code/data/alignments"
#     gff = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
#     config_dict = workflow.process_config("local_config")
#     local = True
#     workflow.run_counts_for_single_genome(gff, bam_folder, config_dict, today, local)
#     file_names = workflow.find_files_in_a_tree(bam_folder, "bam")
#     for file in file_names:
#         count_file = file.split(".bam")[0] + "_counts_not_st.csv"
#         assert os.path.getsize(count_file) != 0
#
#
# def test_run_alignments_for_single_genome(local_mg1655_fastq_folder):
#     genome, fastq_folder, today, config_dict, local = local_mg1655_fastq_folder
#     workflow.run_alignments_for_single_genome(genome, fastq_folder, config_dict, today, local)
#
#     file_names = workflow.find_files_in_a_tree(fastq_folder)  # todo test this function
#
#     for file in file_names:
#         bam_file = file.split(".")[0] + ".bam"
#         assert os.path.isfile(bam_file)
#



if __name__ == "__main__":
    print("Hello!")
