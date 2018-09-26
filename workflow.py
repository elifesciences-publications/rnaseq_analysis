import argparse
import datetime as dt
import os
import pandas as pd
import subprocess
import shlex
import shutil
from modules import generate_PBS_script
from modules import process_fastq
from modules import align
from modules import helpers
from modules import samtools
from modules import counts


def get_args():
    parser = argparse.ArgumentParser("RNASeq Pipeline\n")
    parser.add_argument("-a", "--analysis", help="Analysis options", required=True)
    parser.add_argument("-i", "--input", help="Directory of files", required=True)
    parser.add_argument('-ref', '--reference_genome', help="Reference genome for alignments", required=False)
    parser.add_argument('-gff', '--gff', help="Annotation", required=False)
    parser.add_argument("-local", "--local", help='Run locally or on flux', action='store_true',
                        required=False)
    return parser


def submit_local_job(script):
    multiple_scripts = script.strip().split("\n")
    out = []
    for scr in multiple_scripts:
        cmd = shlex.split(scr)
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output = helpers.to_str(output.stdout.read().strip())
        out.append(output)
    return out


def submit_flux_job(out_dir, sample_id, job_name, script, job_depend=''):
    """

    :param out_dir: where to save .pbs script
    :param sample_id: file_id
    :param job_name: what is being done to that file_id
    :param script:
    :param job_depend:
    :return:
    """
    today = dt.datetime.today().strftime("%Y-%m-%d")
    pbs_name = os.path.join(out_dir, "{}_{}_{}.pbs".format(today, sample_id, job_name))
    generate_PBS_script.generate_PBS_script(pbs_name, script)
    if job_depend:
        output = subprocess.Popen(["qsub", "-W", "depend=afterok:{}".format(job_depend),
                                   pbs_name], stdout=subprocess.PIPE)
    else:
        output = subprocess.Popen(["qsub", pbs_name], stdout=subprocess.PIPE)
    output = output.stdout.read().split(".")[0]
    return helpers.to_str(output)


"""Single File Jobs"""


def run_trim_job(fastq_file, config_dict, local=False, job_depend=''):
    """
    :param fastq_file: file name, if PE forward read, assumes reverse
    is in the same location, with same name except 2 instead of one
    :param config_dict: output of process_config
    :param local: whether job is running locally or on flux
    :param job_depend: if running on flux, whether has to wait for other jobs to finish first
    :return: fastq_file_output path + '' if local, jobid if on flux

    """
    assert '.fastq' in fastq_file
    mode = config_dict["sequencing"]["type"]
    param_dict = config_dict["Trimmomatic"]
    out_dir = os.path.dirname(fastq_file)
    sample_path = helpers.to_str(fastq_file.split(".fastq")[0])
    sample_id = os.path.basename(fastq_file).split(".")[0]
    if mode == "PE":
        # Only handles certain naming formats (see helpers.get_second)
        first = os.path.basename(fastq_file)
        first_out = os.path.join(out_dir, (sample_path + "_trimmed.fastq"))
        second = os.path.join(out_dir, helpers.get_second(first))
        sample_path2 = helpers.to_str(second).split(".fastq")[0]
        second_out = os.path.join(out_dir, (sample_path2 + "_trimmed.fastq"))
        fastq_file = fastq_file + " " + second
        fastq_file_trimmed = first_out + " " + second_out
    else:
        fastq_file_trimmed = sample_path + "_trimmed.fastq"
    script = process_fastq.trimmomatic(fastq_file, fastq_file_trimmed, param_dict)
    if local:
        submit_local_job(script)
        return fastq_file_trimmed, ''
    else:
        jobid = submit_flux_job(out_dir, sample_id, "trimmomatic", script, job_depend)
        return fastq_file_trimmed, jobid


def run_fastqc_job(fastq_file, config_dict, local=False, job_dependency=''):
    param_dict = config_dict["FastQC"]
    sample_id = os.path.basename(fastq_file).split(".")[0]
    out_dir = os.path.dirname(fastq_file)
    fastqc_out_dir = os.path.join(out_dir, "{}_FastQC_results".format(sample_id))
    subprocess.call(["mkdir", "-p", fastqc_out_dir])
    script = process_fastq.fastqc(fastq_file, fastqc_out_dir, param_dict)
    if local:
        submit_local_job(script)
        return fastqc_out_dir, ''
    else:
        jobid = submit_flux_job(out_dir, sample_id,  "FastQC", script, job_dependency)
        return fastqc_out_dir, jobid


def run_multiqc_job(fastqc_dir, config_dict, local, job_depend=''):
    # todo test
    param_dict = config_dict["MultiQC"]
    report_name = dt.datetime.today().strftime("%Y-%m-%d") + "multiqc_report"
    script = process_fastq.multiqc(fastqc_dir, fastqc_dir, param_dict)
    if local:
        submit_local_job(script)
        return fastqc_dir, ""
    else:
        jobid = submit_flux_job(fastqc_dir, report_name,  "Mqc", script, job_depend)
        return fastqc_dir, jobid


def run_build_index_job(ref_genome, config_dict, local=False, job_depend=''):
    param_dict = config_dict["Bowtie"]
    ref_path = helpers.to_str(ref_genome.split(".")[0])
    ref_id = os.path.basename(ref_genome).split(".")[0]
    index_path = ref_path + '_index'
    script = align.build_bowtie_index(ref_genome, index_path, param_dict)
    if local:
        submit_local_job(script)
        return index_path, ''
    else:
        out_dir = os.path.dirname(ref_genome)
        jobid = submit_flux_job(out_dir, ref_id, "bowtie_index", script, job_depend)
        return index_path, jobid
    # todo test index job on flux


def run_alignment_job(fastq_file, index_path, config_dict, local, job_depend=''):
    out_dir = os.path.dirname(fastq_file)
    sample_id = helpers.to_str(os.path.basename(fastq_file).split(".")[0])
    param_dict = config_dict["Bowtie"]
    seq_type = config_dict["sequencing"]["type"]
    if seq_type == "PE":
        second = helpers.get_second(fastq_file)
        fastq_file = fastq_file + " " + second
    sam_file = os.path.join(out_dir, sample_id + ".sam")
    script = align.bowtie_align(fastq_file, sam_file, index_path, param_dict)
    if local:
        submit_local_job(script)
        return sam_file, ''
    else:
        jobid = submit_flux_job(out_dir, sample_id, "bowtie_align", script, job_depend)
        return sam_file, jobid


def run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local, job_dependency=''):
    param_dict = config_dict["Samtools"]
    out_dir = os.path.dirname(sam_file)
    suffix = sam_file.split(".")[0]
    bam_file = suffix + ".bam"
    sorted_bam_file = bam_file.split(".bam")[0]+"_sorted.bam"
    script = samtools.sam2bam(sam_file, bam_file, param_dict)
    if local:
        submit_local_job(script)
        return sorted_bam_file, ''
    else:
        suffix = os.path.basename(sam_file).split(".")[0]
        jobid = submit_flux_job(out_dir, suffix, "sam2bam", script, job_dependency)
        return sorted_bam_file, jobid


def run_bam_stats(bam_file):
    total, mapped, pnt_mapped = helpers.get_bam_stats(bam_file)
    return total, mapped, pnt_mapped


def run_count_job_bedtools(gff, bam, config_dict, today, local, job_dependency=""):
    param_dict = config_dict["bedtools"]
    output_directory = os.path.dirname(bam)
    prefix = os.path.basename(bam).split(".")[0]
    strand = True if config_dict["bedtools"]["strand"] == "-s" else False
    suffix = "st" if strand else "not_st"
    count_file = bam.split(".bam")[0] + "_counts_{}.csv".format(suffix)
    if local:
        script = counts.count_with_bedtools(gff, bam, param_dict)
        script = script.split(">")[0]
        result = submit_local_job(script)[0]
        with open(count_file, "w") as fo:
            fo.write(result)
        return count_file, ''
    else:
        script = counts.count_with_bedtools(gff, bam, param_dict)
        jobid = submit_flux_job(output_directory, prefix, "Count", script, job_dependency)
        return count_file, jobid


def run_edit_count_job_bedtools(count_file, config_dict, today, local, job_dependency=""):
    param_dict = config_dict["bedtools"]
    count_file_edited = count_file.split(".csv")[0] + "_edited.csv"
    out_dir = os.path.dirname(count_file)
    prefix = os.path.basename(count_file).split(".")[0]
    if local:
        helpers.process_bedtools_count_output(count_file, count_file_edited, param_dict)
        return count_file_edited, ""
    else:
        # todo make this less hacky
        # todo test
        config_file = "/Users/annasintsova/git_repos/code/local_config"
        script = "python modules/counts.py {} {} {}".format(count_file, count_file_edited, config_file)
        jobid = submit_flux_job(out_dir, prefix, "count_edit", script, job_dependency)
        return count_file_edited, jobid



def run_count_job_htseq_count(gff, bam, config_dict, today, local, job_dependency=""):
    param_dict = config_dict["HTSeq"]
    script = counts.count_reads(bam, gff, param_dict)


def clean_up_after_run(fastq_dir, out_dir):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    bam_dir = os.path.join(out_dir, "bams")
    counts_dir = os.path.join(out_dir, "counts")
    os.mkdir(bam_dir)
    os.mkdir(counts_dir)
    files = [f for f in os.listdir(fastq_dir)]
    for f in files:
        if ".bam" in f:
            shutil.move(os.path.join(fastq_dir, f), os.path.join(bam_dir, f))
        elif "counts" in f:
            shutil.move(os.path.join(fastq_dir, f), os.path.join(counts_dir, f))
        elif ".sam" in f:
            os.remove(os.path.join(fastq_dir, f))


"""Workflows"""


def workflow_test(analysis, input_folder, output_folder):
    return analysis, input_folder, output_folder


def workflow_qc(fastq_dir, config_dict, local=False):
    files = helpers.find_files_in_a_tree(fastq_dir, file_type='fastq')
    for file in files:
        run_fastqc_job(file, config_dict, local)
    return "Fastqc jobs submitted!"


def workflow_mqc(fastqc_dir, config_dict, local):
    run_multiqc_job(fastqc_dir, config_dict, local)
    return "Multiqc jobs submitted!"

# todo make sure can handle .gz files


def workflow_trim_and_qc(fastq_dir, config_dict, local=False):

    files = helpers.find_files_in_a_tree(fastq_dir, file_type='.fastq')
    for file in files:
        # Run trim job
        trimmed_fastq_file, trim_jobid = run_trim_job(file, config_dict, local)
        # Run fastqc job
        run_fastqc_job(trimmed_fastq_file, config_dict, local, trim_jobid)
    return "Jobs for trim and qc submitted"  # todo test workflow on flux


def workflow_align(fastq_dir, ref, gff, config_dict, today, local):

    if local:
        # 1. Build index
        print("Building index...")
        bt2, _ = run_build_index_job(ref, today, config_dict, local)
        print("Index complete, index name: {}".format(bt2))
        # 2. Find fastq files
        print("Looking for fastq files")
        fastq_files = helpers.find_files_in_a_tree(fastq_dir, file_type='fastq')
        print("Found {} fastq files".format(len(fastq_files)))
        # 3. Iterate over them and align
        for file in fastq_files:
            print("Aligning {}".format(file))
            sam_file, _ = run_alignment_job(file, bt2, config_dict, today, local)
            print("Sorting {}".format(sam_file))
            sorted_bam, _ = run_sam_to_bam_conversion_and_sorting(sam_file, config_dict, today, local)
            print("Counting {}".format(sorted_bam))
            counts_file, _ = run_count_job_bedtools(gff, sorted_bam, config_dict, today, local)
            print("Counting complete, count file: {}".format(counts_file))
            print("Editing count file")
            run_edit_count_job_bedtools(counts_file, config_dict, today, local)

    else:
        job_ids = {}
        # 1. Build index
        bt2, index_jobid = run_build_index_job(ref, today, config_dict, local)
        # 2. Find fastq files
        fastq_files = helpers.find_files_in_a_tree(fastq_dir, file_type='fastq')
        # 3. Iterate over them and align
        for file in fastq_files:
            sam_file, samfile_jobid = run_alignment_job(file, bt2, config_dict,
                                                        today, local, index_jobid)
            sorted_bam, sam2bam_jobid = run_sam_to_bam_conversion_and_sorting(sam_file, config_dict,
                                                                              today, local,
                                                                              samfile_jobid)
            counts_file, counts_jobid = run_count_job_bedtools(gff, sorted_bam, config_dict, today,
                                                               local, sam2bam_jobid)
            count_file_edited, ce_jobid = run_edit_count_job_bedtools(counts_file, config_dict, today,
                                                                      local, job_dependency=counts_jobid)

            job_ids[file] = [samfile_jobid, sam2bam_jobid, counts_jobid, ce_jobid]
        return job_ids


def workflow_bam_stats(bam_dir, today, local, job_dependency=""):
    if local:
        bam_files = helpers.find_files_in_a_tree(bam_dir, file_type="bam")
        stats = []
        labels = ["Name", "Total", "Mapped", "% Mapped"]
        for bm in bam_files:
            sample_name = os.path.basename(bm).split(".")[0]
            total, mapped, pcnt_mapped = run_bam_stats(bm)
            stats.append((sample_name, total, mapped, pcnt_mapped))
        df = pd.DataFrame.from_records(stats, index="Name", columns=labels)
        filename = os.path.join(bam_dir, today+"_alignment_stats.csv")
        df.to_csv(filename)
    else:
        script = "python workflow.py -a bam-stats -i {} -local".format(bam_dir)
        jobid = submit_flux_job(bam_dir, "alignment_stats",  "alignment_stats", script, job_dependency)
        return jobid



# def workflow_align(genome, fastq_folder, config_dict, today, local):
#
#     run_alignments_for_single_genome(genome, fastq_folder, config_dict, today, local)
#     return "Simple Align Workflow!"


def workflow_count(gff, bam_folder, config_dict, today, local):
    bams = helpers.find_files_in_a_tree(bam_folder, file_type="bam")
    if local:
        for bam in bams:
            count_file, _ = run_count_job_bedtools(gff, bam, config_dict, today, local)
            count_file_edited, _ = run_edit_count_job_bedtools(count_file, config_dict, today, local)
    else:
        for bam in bams:
            count_file, jobid = run_count_job_bedtools(gff, bam, config_dict, today, local)
            count_file_edited, _ = run_edit_count_job_bedtools(count_file, config_dict, today, local,
                                                               job_dependency=jobid)


# def run_count_job_htseq(gff, bam, config_dict, local, job_dependency=""):
#


# def run_alignments_for_multiple_genomes(genome_read_pairs, today, config_dict): #list of tuples
#     bowtie_bin = config_dict["Bowtie"]["bin"]
#     samtools_bin = config_dict["Samtools"]["bin"]
#     for genome_read_pair in genome_read_pairs:
#         genome = genome_read_pair[0]
#         fastq_file = genome_read_pair[1]
#         output_directory = genome_read_pair[2]
#         #run build index
#         bt2_base, index_jobid = run_build_index_job(genome, output_directory,
#                             today, bowtie_bin)
#         #run alignment job
#         sam_file, align_jobid = run_alignment_job(fastq_file, output_directory,
#                           bt2_base, bowtie_bin, today,
#                           job_dependency=index_jobid)
#         #run sam job
#         # bam_file, samtools_jobid = run_sam_to_bam_conversion_and_sorting(sam_file,
#         #                                                                  output_directory,
#         #                                                                  today,
#         #                                                                  samtools_bin,
#         #                                                                  job_dependency=align_jobid)

#

def flow_control():
    today = dt.datetime.now().strftime("%Y-%m-%d")
    args = get_args().parse_args()

    if args.local:
        config_dict = helpers.process_config("local_config")
    else:
        config_dict = helpers.process_config("config")
    if args.analysis == 'test':
        print(workflow_test(args.analysis, args.input, args.out_dir))
    elif args.analysis == 'qc':
        print(workflow_qc(args.input, config_dict, args.local))
    elif args.analysis == 'mqc':
        print(workflow_mqc(args.input, config_dict, args.local))
    elif args.analysis == 'trim':
        print(workflow_trim_and_qc(args.input, config_dict, args.local))
    elif args.analysis == 'align':  # tested locally
        print(workflow_align(args.input, args.reference_genome, args.gff, config_dict, today, args.local))
    elif args.analysis == 'bam-stats':
        print(workflow_bam_stats(args.input, today, args.local))
    elif args.analysis == 'count':  # todo test
        assert args.gff
        workflow_count(args.gff, args.input, config_dict, today, args.local)


if __name__ == "__main__":
    flow_control()
