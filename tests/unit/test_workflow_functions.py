import datetime as dt
import os
import pytest
import subprocess
import sys
sys.path.append('.')
import workflow
from modules.helpers import to_str


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
    today = dt.datetime.today().strftime("%Y-%m-%d")
    job_name = "test_job"
    script = "mkdir TEST1"


    output = workflow.submit_flux_job(output_directory, suffix, today,
                             job_name, script)
    #assert type(output) in [str, unicode] # For python 2 (Flux)
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




if __name__ == "__main__":
    print("Hello!")
