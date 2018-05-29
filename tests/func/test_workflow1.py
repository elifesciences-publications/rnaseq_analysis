import sys
sys.path.append('.')
import subprocess
import shlex
import workflow
from modules.helpers import to_str

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

def test_workflow1():

    # User enters location of fastq files, and output directory,
    # Run default trimmomatic
    # Run default fastqc/multiqc, output in the output directory

    analysis = 'workflow1'
    input_folder = "/Users/annasintsova/git_repos/code/data/reads"
    output_folder = "/Users/annasintsova/git_repos/code/tests/test_data"
    cmd_str = "python workflow.py -a {} -i {} -o {}".format(analysis,
                                                            input_folder, output_folder)
    cmd = shlex.split(cmd_str)
    output1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output1 = to_str(output1.stdout.read())

    # What am I asserting?
    # PBS script creation first
    # Eventually the right files in right locations? The correctness of files will be unit tests?

    assert(output1 == "Jobs for workflow1 submitted\n")



if __name__ == "__main__":
    test_can_call_workflow_with_arguments()