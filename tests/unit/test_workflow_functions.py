import datetime as dt
import os
import pytest
import sys
sys.path.append('.')
import workflow
from modules.helpers import to_str


# today = dt.datetime.today().strftime("%Y-%m-%d")

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


@pytest.mark.skip(reason="only on FLUX")
def test_submit_flux_job_one_job():
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data/"
    suffix = "test_script1"
    today = dt.datetime.today().strftime("%Y-%m-%d")  # todo refactor today variable
    job_name = "test_job"
    script = "mkdir TEST1"
    output = workflow.submit_flux_job(output_directory, suffix, today,
                                      job_name, script)
    assert type(int(output)) == int


@pytest.mark.skip(reason="only on FLUX")
def test_submit_flux_job_two_jobs():
    # Submit first job:
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data/"
    suffix = "test_script1"
    today = dt.datetime.today().strftime("%Y-%m-%d")
    job_name = "test_job"
    script = "mkdir TEST1"

    output1 = workflow.submit_flux_job(output_directory, suffix, today,
                                       job_name, script)

    # Submit second job:
    suffix2 = 'test_script2'
    script2 = "cd TEST1\nmkdir TEST2"

    output2 = workflow.submit_flux_job(output_directory, suffix2, today,
                                       job_name, script2, job_dependency=output1)
    assert type(int(output2)) == int


# Both of these will test local system

def test_run_trim_job(local_fastq_hm86_ur_tmpdir_out):

    # GIVEN:
    filename, outdir, td, config_dict, loc = local_fastq_hm86_ur_tmpdir_out
    # WHEN:
    output_file_name = workflow.run_trim_job(filename, outdir, td, config_dict, loc)[0]
    # THEN:
    assert os.path.isfile(output_file_name)
    assert os.path.getsize(output_file_name) != 0


def test_run_fastqc_job(local_fastq_hm86_ur_tmpdir_out):
    # GIVEN:
    filename, outdir, td, config_dict, loc = local_fastq_hm86_ur_tmpdir_out
    # WHEN:
    actual_out_dir = workflow.run_fastqc_job(filename, outdir, td, config_dict, loc)[0]
    # THEN:
    expected_filename = os.path.join(actual_out_dir,
                                     to_str(os.path.basename(filename).split(".fastq")[0]) + "_fastqc.html")

    assert os.path.join(outdir, "FastQC_results") == actual_out_dir
    assert len(os.listdir(actual_out_dir)) != 0
    assert os.path.isfile(expected_filename)

# todo make tests clean up after themselves

# this will test both on flux, plus if job dependency works properly


@pytest.mark.skip(reason="only on FLUX")
def test_run_fastqc_after_run_trim_job():
    fastq_file_input = "/scratch/hmobley_fluxod/annasint/code/data/reads/HM86_UR_copy.fastq.gz"
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data"
    today = dt.datetime.today().strftime("%Y-%m-%d")
    config_dict = workflow.process_config(config_file="config")
    local = False
    output_file_name, jobid = workflow.run_trim_job(fastq_file_input, output_directory,
                                                    today, config_dict, local)
    assert type(int(jobid)) == int  # only way can figure out if the job has been submitted
    actual_out_dir, jobid2 = workflow.run_fastqc_job(output_file_name, output_directory, today,
                                                     config_dict, local, job_dependency=jobid)
    assert type(int(jobid2)) == int  # todo come up with a better assert statement


def test_run_build_index_job(local_mg1655_index):
    reference_genome, fastq_file, today, config_dict, local = local_mg1655_index
    bt2, _ = workflow.run_build_index_job(reference_genome,
                                          today, config_dict,
                                          local)
    assert bt2 == "/Users/annasintsova/git_repos/code/data/ref/MG1655_index"


def test_run_align_job_local(local_mg1655_index):
    reference_genome, fastq_file, today, config_dict, local = local_mg1655_index # todo make sure sam file exists
    bt2, _ = workflow.run_build_index_job(reference_genome,
                                          today, config_dict,
                                          local)
    sam_file, _ = workflow.run_alignment_job(fastq_file,
                                          bt2, config_dict, today,
                                          local)
    assert os.path.isfile(sam_file)

def test_run_sam_to_bam_conversion_and_sorting_local(local_mg1655_index):
    sam_file = "/Users/annasintsova/git_repos/code/data/reads/SRR1051511_trimmed.sam"
    _, _, today, config_dict, local = local_mg1655_index

    bam_file, _ = workflow.run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local)
    expected_bam_sorted = "/Users/annasintsova/git_repos/code/data/reads/SRR1051511_trimmed_sorted.bam"
    expected_index_bam = "/Users/annasintsova/git_repos/code/data/reads/SRR1051511_trimmed_sorted.bam.bai"
    expected_flagstat_file = "/Users/annasintsova/git_repos/code/data/reads/SRR1051511_trimmed_flagstat.txt"
    assert os.path.isfile(bam_file)
    assert os.path.isfile(expected_bam_sorted)
    assert os.path.isfile(expected_index_bam)
    #assert os.path.isfile(expected_flagstat_file)

def test_run_alignments_for_single_genome(local_mg1655_fastq_folder):
    genome, fastq_folder, today, config_dict, local = local_mg1655_fastq_folder
    workflow.run_alignments_for_single_genome(genome, fastq_folder, config_dict, today, local)

    file_names = workflow.find_files_in_a_tree(fastq_folder)  # todo test this function

    for file in file_names:
        bam_file = file.split(".")[0] + ".bam"
        assert os.path.isfile(bam_file)

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


if __name__ == "__main__":
    print("Hello!")
