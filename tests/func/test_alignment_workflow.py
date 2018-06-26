import sys
sys.path.append('.')
from modules.helpers import to_str
import os
import subprocess
import shlex

# User enters location of fastq files, and output directory,
# Run default trimmomatic
# Run default fastqc/multiqc, output in the output directory


def test_alginment_workflow_locally():

    analysis = 'align'
    input_folder = "/Users/annasintsova/git_repos/code/data/reads"
    local = '-local'
    ref = "/Users/annasintsova/git_repos/code/data/ref/MG1655.fna"
    cmd_str = "python workflow.py -a {} -i {} -ref {} {}".format(analysis,
                                                                 input_folder,
                                                                 ref,
                                                              local)
    cmd = shlex.split(cmd_str)

    subprocess.call(cmd)

    expected_files = ["/Users/annasintsova/git_repos/code/data/reads/HM_fastq/SRR1051490.fastq",
                      "/Users/annasintsova/git_repos/code/data/reads/SRR_fastq/SRR1051514.fastq"]

    for fi in expected_files:
        suffix = fi.split(".fastq")[0]
        assert os.path.isfile(suffix + ".bam")
        suffix2 = suffix + "_sorted.bam"
        assert os.path.isfile(suffix2)
        assert os.path.isfile(os.path.join(input_folder,  'bam_alignment_stats.csv') )

#
# if __name__ == "__main__":
#     test_alginment_workflow_locally()