from modules.helpers import to_str
import os
import subprocess
import shlex

# A user can call workflow.py and pass in some parameters


def test_can_call_workflow_with_arguments():

    analysis = 'test'
    input_folder = "/Users/annasintsova/git_repos/code/data/reads"
    output_folder = "/Users/annasintsova/git_repos/code/tests/test_data"
    cmd_str = "python workflow.py -a {} -i {} -o {}".format(analysis,
                                                            input_folder, output_folder)
    cmd = shlex.split(cmd_str)

    output1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output1 = to_str(output1.stdout.read())

    assert(analysis in output1)
    assert (input_folder in output1)
    assert(output_folder in output1)


def test_workflow1_locally():

    # User enters location of fastq files, and output directory,
    # Run default trimmomatic
    # Run default fastqc/multiqc, output in the output directory

    analysis = 'workflow1'
    input_folder = "/Users/annasintsova/git_repos/code/data/reads"
    output_folder = "/Users/annasintsova/git_repos/code/tests/test_data"
    local = '-local'
    cmd_str = "python workflow.py -a {} -i {} -o {} {}".format(analysis,
                                                               input_folder,
                                                               output_folder, local)
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    expected_files = [fi for fi in os.listdir(input_folder)]

    for fi in expected_files:

        suffix = fi.split(".fastq")[0]
        out_dir = os.path.join(output_folder, suffix+"_output")
        assert os.path.isfile(os.path.join(out_dir, suffix+"_trimmed.fastq"))
        suffix2 = fi.split(".fastq")[0] + "_trimmed" + "_fastqc.html"
        assert os.path.isfile(os.path.join(out_dir, "FastQC", suffix2))


if __name__ == "__main__":
    test_can_call_workflow_with_arguments()
