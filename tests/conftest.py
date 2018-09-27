import datetime as dt
import os
import pytest
import sys
import shutil
sys.path.append('.')
import workflow
from modules import helpers

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
    config_dict = helpers.process_config(config_file="local_config")
    local = True
    yield (str(f_file), str(r_file), config_dict, local)

@pytest.fixture()
def local_sam(tmpdir, day):
    sam_file = "/Users/annasintsova/git_repos/code/data/alignments/SRR1051490.sam"
    s_file = tmpdir.mkdir("alignments").join("align.sam")
    shutil.copy(sam_file, str(s_file))
    config_dict = helpers.process_config(config_file="local_config")
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
    config_dict = helpers.process_config(config_file="local_config")
    local = True
    today = day
    yield (str(b_file), str(g_file), today, config_dict, local)


@pytest.fixture()
def local_fastq_dir_ref(tmpdir, day):

    data_dir = "/Users/annasintsova/git_repos/code/data/reads"
    files = helpers.find_files_in_a_tree(data_dir, "fastq")
    for f in files:
        fname = os.path.basename(f)
        new_name = str(tmpdir.join(fname))
        shutil.copy(f, new_name)
    ref = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
    ref_name = os.path.basename(ref)
    test_ref = str(tmpdir.join(ref_name))
    shutil.copy(ref, test_ref)
    return str(tmpdir), test_ref, day


"""FLUX"""

@pytest.fixture()
def flux_fastq_ref(tmpdir):
    fastq_file = "/scratch/hmobley_fluxod/annasint/code/data/reads/SRR1051490.fastq"
    reference_genome = "/scratch/hmobley_fluxod/annasint/code/data/ref/MG1655.fna"
    #f_file = tmpdir.join("test.fastq")
    #r_file = tmpdir.mkdir("ref").join("genome.fna")
    test_dir = "/scratch/hmobley_fluxod/annasint/code/test_data"
    os.mkdir(test_dir)
    f_file = os.path.join(test_dir, "test.fastq")
    r_file = os.path.join(test_dir, "genome.fna")
    shutil.copy(fastq_file, f_file)
    shutil.copy(reference_genome, r_file)
    config_dict = helpers.process_config(config_file="config")
    local = False
    yield (f_file, r_file, config_dict, local)


@pytest.fixture()
def flux_fastq_dir():
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/reads/"
    f_dir = "/scratch/hmobley_fluxod/annasint/code/test_data"
    shutil.copytree(data_dir, f_dir)
    config_dict = helpers.process_config(config_file="config")
    local = False
    yield (f_dir, config_dict, local)


@pytest.fixture()
def flux_mqc_dir():
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/mqc_dir/"
    mqc_dir = "/scratch/hmobley_fluxod/annasint/code/test_data"
    shutil.copytree(data_dir, mqc_dir)
    config_dict = helpers.process_config(config_file="config")
    local = False
    yield (mqc_dir, config_dict, local)


@pytest.fixture()
def flux_sam_bam(tmpdir, day):
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/"

    gff = os.path.join(data_dir, "ref/MG1655.gff")
    gff_name = "MG1655.gff"
    test_gff = str(tmpdir.join(gff_name))
    shutil.copy(gff, test_gff)

    sam = os.path.join(data_dir, "alignments/SRR1051490.sam")
    sam_name = "SRR1051490.sam"
    test_sam = str(tmpdir.join(sam_name))
    shutil.copy(sam, test_sam)

    bam = os.path.join(data_dir, "alignments/SRR1051490_sorted.bam")
    bam_name = "SRR1051490_sorted.bam"
    test_bam = str(tmpdir.join(bam_name))
    shutil.copy(bam, test_bam)

    config_dict = helpers.process_config(config_file="config")
    local = False
    return test_sam, test_bam, test_gff, config_dict, local

@pytest.fixture()
def flux_bam_keep_output():
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/"
    test_dir = "/scratch/hmobley_fluxod/annasint/code/test_data"
    os.mkdir(test_dir)

    gff = os.path.join(data_dir, "ref/MG1655.gff")
    test_gff = os.path.join(test_dir, "MG1655.gff")
    shutil.copy(gff, test_gff)

    bam = os.path.join(data_dir, "alignments/SRR1051490_sorted.bam")
    test_bam = os.path.join(test_dir, "SRR1051490_sorted.bam")
    shutil.copy(bam, test_bam)
    bai = os.path.join(data_dir, "alignments/SRR1051490_sorted.bam.bai")
    test_bai = os.path.join(test_dir, "SRR1051490_sorted.bam.bai")
    shutil.copy(bai, test_bai)
    config_dict = helpers.process_config(config_file="config")
    local = False
    return test_bam, test_gff, config_dict, local


@pytest.fixture()
def flux_fastq_dir_ref(tmpdir, day):
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/reads"
    files = helpers.find_files_in_a_tree(data_dir, "fastq")
    for f in files:
        fname = os.path.basename(f)
        new_name = str(tmpdir.join(fname))
        shutil.copy(f, new_name)
    ref = "/scratch/hmobley_fluxod/annasint/code/data/ref/MG1655.fna"
    ref_name = os.path.basename(ref)
    test_ref = str(tmpdir.join(ref_name))
    shutil.copy(ref, test_ref)
    return str(tmpdir), test_ref, day


@pytest.fixture()
def count_keep_output():
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data/"
    test_dir = "/scratch/hmobley_fluxod/annasint/code/test_data"
    os.mkdir(test_dir)

    cnts = os.path.join(data_dir, "counts/SRR1051490_sorted_counts_st.csv")
    test_cnts = os.path.join(test_dir, "SRR1051490_sorted_counts_st.csv")
    shutil.copy(cnts, test_cnts)
    config_dict = helpers.process_config(config_file="config")
    local = False
    return test_cnts, config_dict, local


@pytest.fixture()
def flux_fastq_dir_ref_PE(day): # todo refactor
    data_dir = "/scratch/hmobley_fluxod/annasint/code/data2/"
    test_dir = "/scratch/hmobley_fluxod/annasint/code/test_data/"
    shutil.copytree(data_dir, test_dir)
    ref_genome = os.path.join(test_dir, 'ref/MG1655.fna')
    gff = os.path.join(test_dir, 'ref/MG1655.gff')
    fastq_file_dir = os.path.join(test_dir, "reads")
    today = day
    config_dict = helpers.process_config(config_file="config")
    local = False
    return fastq_file_dir, ref_genome, gff, today, config_dict, local





# @pytest.fixture()
# def remove_test_data():
#     yield
#     os.remove()
# todo finish fixture that would remove all data from test_data directory in the end