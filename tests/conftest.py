import datetime as dt
import os
import pytest
import sys
sys.path.append('.')
import workflow

@pytest.fixture()
def day():
    return dt.datetime.today().strftime("%Y-%m-%d")


@pytest.fixture()
def local_fastq_hm86_ur_tmpdir_out(tmpdir, day):
    fastq_file = "/Users/annasintsova/git_repos/code/data/reads/HM86_UR_copy.fastq.gz"
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    yield (fastq_file, str(tmpdir), today, config_dict, local)

@pytest.fixture()
def local_mg1655_index(day):
    reference_genome = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
    fastq_file = "/Users/annasintsova/git_repos/code/"\
                 "data/reads/SRR1051511_trimmed.fastq"
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    return reference_genome, fastq_file, today, config_dict, local


@pytest.fixture()
def local_mg1655_fastq_folder(day):
    reference_genome = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
    fastq_folder = "/Users/annasintsova/git_repos/code/data/reads/"
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    return reference_genome, fastq_folder, today, config_dict, local


# @pytest.fixture()
# def remove_test_data():
#     yield
#     os.remove()
# todo finish fixture that would remove all data from test_data directory in the end