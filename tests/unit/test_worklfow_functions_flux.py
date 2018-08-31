import datetime as dt
import os
import pytest
import sys
import shutil
import workflow
from modules.helpers import to_str

sys.path.append('.')


def test_submit_flux_job_one_job():
    output_directory = "/scratch/hmobley_fluxod/annasint/code/tests/test_data/"
    suffix = "test_script1"
    today = dt.datetime.today().strftime("%Y-%m-%d")  # todo refactor today variable
    job_name = "test_job"
    script = "mkdir TEST1"
    output = workflow.submit_flux_job(output_directory, suffix, today,
                                      job_name, script)
    assert type(int(output)) == int


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

# todo make tests clean up after themselves

# this will test both on flux, plus if job dependency works properly


def test_run_fastqc_job(flux_fastq_ref):
    fastq_file_input, _, today, config_dict, local = flux_fastq_ref
    output_file_name, jobid = workflow.run_fastqc_job(fastq_file_input, today, config_dict, local)
    assert type(int(jobid)) == int  # only way can figure out if the job has been submitted


def test_run_trim_job(flux_fastq_ref):
    fastq_file_input, _,  today, config_dict, local = flux_fastq_ref
    output_file_name, jobid = workflow.run_trim_job(fastq_file_input, today, config_dict, local)
    assert type(int(jobid)) == int


def test_run_index_job(flux_fastq_ref):
    fastq_file_input, ref, today, config_dict, local = flux_fastq_ref
    bt2, jobid = workflow.run_build_index_job(ref, today, config_dict, local)
    assert type(int(jobid)) == int


def test_run_fastqc_after_run_trim_job(flux_fastq_ref):
    fastq_file_input, _, today, config_dict, local = flux_fastq_ref
    output_file_name, jobid = workflow.run_trim_job(fastq_file_input, today, config_dict, local)
    assert type(int(jobid)) == int  # only way can figure out if the job has been submitted
    actual_out_dir, jobid2 = workflow.run_fastqc_job(output_file_name, today, config_dict, local, job_dependency=jobid)
    assert type(int(jobid2)) == int  # todo come up with a better assert statement


def test_run_align_job(flux_fastq_ref):
    fastq_file, ref, today, config_dict, local = flux_fastq_ref
    bt2, jobid = workflow.run_build_index_job(ref, today, config_dict, local)
    assert type(int(jobid)) == int
    sam_file, jobid2 = workflow.run_alignment_job(fastq_file, bt2, config_dict, today, local, jobid)
    assert type(int(jobid2)) == int

def test_run_samtobam_job(flux_sam_bam):
    sam_file, _, _,  today, config_dict, local = flux_sam_bam
    bam_file, jobid = workflow.run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local)
    assert type(int(jobid)) == int


def test_run_count_job(flux_sam_bam):
    _,  bam_file, gff_file, today, config_dict, local = flux_sam_bam
    counts_file, jobid = workflow.run_count_job_bedtools(gff_file, bam_file, config_dict, today, local)
    assert type(int(jobid)) == int


def test_worflow2a(flux_fastq_dir_ref):
    fastq_file_dir, ref_genome, gff, today, config_dict, local = flux_fastq_dir_ref
    jobids_dict = workflow.workflow2a(ref_genome, fastq_file_dir, gff, config_dict, today, local)
    fastq_files = workflow.find_files_in_a_tree(fastq_file_dir)
    for f in fastq_files:
        for v in jobids_dict[f]:
            assert type(int(v)) == int


def test_workflow2a_PE(flux_fastq_dir_ref_PE):
    fastq_file_dir, ref_genome, gff, today, config_dict, local = flux_fastq_dir_ref_PE
    jobids_dict = workflow.workflow2a(ref_genome, fastq_file_dir, gff, config_dict, today, local)
    for key, val in jobids_dict.items():
        print(key)
        for v in val:
            assert type(int(v)) == int








