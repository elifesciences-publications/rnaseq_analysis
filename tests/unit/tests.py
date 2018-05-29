import datetime as dt
import os
import pytest
import subprocess
import sys
sys.path.append('.')
import workflow
from modules import process_fastq
from modules import generate_PBS_script as pbs
from modules.helpers import to_str


def test_process_config():
    output = workflow.process_config(config_file="tests/test_data/config")
    desired_output = {'bin_path': {'binbase': '/home/annasint/'}}
    assert output == desired_output



def test_submit_job(pbs_name, job_dependency=''):

    if job_dependency == '1':
        output = '2'
    else:
        output = subprocess.Popen(["cat", pbs_name],stdout=subprocess.PIPE)
        output = output.stdout.read().strip()
    return to_str(output)


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

#############################################################################################################
# MODULES TESTS


def test_modules_trimmomatic():
    input_file = "/Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq"
    output_directory = "/Users/annasintsova/git_repos/code/tests/test_data"
    suffix = os.path.basename(input_file).split(".fastq")[0]
    fastq_file_output = os.path.join(output_directory, suffix + "_trimmed.fastq")
    trimmomatic_bin = "/Users/annasintsova/tools/Trimmomatic-0.36/trimmomatic-0.36.jar"
    trimmomatic_adapters = "/Users/annasintsova/tools/Trimmomatic-0.36/adapters/TruSeq3-SE.fa"

    script = process_fastq.Trimmomatic(input_file, fastq_file_output,
                                       trimmomatic_bin, trimmomatic_adapters)
    assert script == "java -jar /Users/annasintsova/tools/Trimmomatic-0.36/trimmomatic-0.36.jar " \
                     "SE /Users/annasintsova/git_repos/code/data/reads/UTI24_control.fastq " \
                     "/Users/annasintsova/git_repos/code/tests/test_data/UTI24_control_trimmed.fastq " \
                     "ILLUMINACLIP:/Users/annasintsova/tools/" \
                     "Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10:8:true" \
                     " SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0\n"



def test_generate_pbs():

    script_name = "/Users/annasintsova/git_repos/code/tests/test_data/test_script.pbs"
    script = "mkdir TEST1"

    pbs.generate_PBS_script(script_name, script)

    with open(script_name, "r") as fo:
        text = fo.read()
        print(text)
        assert script in text
        assert "#PBS -l nodes=1:ppn=4,pmem=4gb,walltime=12:00:00" in text


# PBS -l nodes=1:ppn=4,pmem={},walltime={}


if __name__ == "__main__":
    test_generate_pbs()
