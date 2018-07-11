import datetime as dt
import os
import pytest
import sys
import shutil
sys.path.append('.')
import workflow

@pytest.fixture()
def day():
    return dt.datetime.today().strftime("%Y-%m-%d")


@pytest.fixture()
def local_fastq_ref(tmpdir, day):
    fastq_file = "/Users/annasintsova/git_repos/code/data/reads/SRR1051490.fastq"
    reference_genome = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
    f_file = tmpdir.join("test.fastq")
    r_file = tmpdir.mkdir("ref").join("genome.fna")
    shutil.copy(fastq_file, str(f_file))
    shutil.copy(reference_genome, str(r_file))
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    yield (str(f_file), str(r_file), today, config_dict, local)

@pytest.fixture()
def local_sam(tmpdir, day):
    sam_file = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490.sam"
    s_file = tmpdir.mkdir("alignments").join("align.sam")
    shutil.copy(sam_file, str(s_file))
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    yield (str(s_file), today, config_dict, local)

@pytest.fixture()
def local_bam(tmpdir, day):
    gff_file = "/Users/annasintsova/git_repos/code/data/ref/MG1655.gff"
    bam_file = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490_sorted.bam"
    bai_file = bam_file + ".bai"
    b_file = tmpdir.mkdir("alignments").join("align.bam")
    bi_file = tmpdir.join("alignments/align.bam.bai")
    g_file = tmpdir.mkdir("ref").join("annotation.gff")
    shutil.copy(bam_file, str(b_file))
    shutil.copy(bai_file, str(bi_file))
    shutil.copy(gff_file, str(g_file))
    config_dict = workflow.process_config(config_file="local_config")
    local = True
    today = day
    yield (str(b_file), str(g_file), today, config_dict, local)


# @pytest.fixture()
# def local_mg1655_fna(tmpdir, day):
#     reference_genome = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
#     out_file = tmpdir.mkdir("ref").join("test.fastq")
#     config_dict = workflow.process_config(config_file="local_config")
#     local = True
#     today = day
#     return reference_genome, fastq_file, today, config_dict, local


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