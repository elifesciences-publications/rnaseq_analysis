import sys
sys.path.append('.')
from modules.helpers import to_str
import os
import subprocess
import shlex
import workflow

def test_alginment_workflow_locally(local_fastq_dir_ref):
    fastq_dir, ref, today = local_fastq_dir_ref
    local = '-local'
    cmd_str = "python workflow.py -a align -i {} -ref {} {}".format(fastq_dir, ref, local)
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)

    fastq_files = workflow.find_files_in_a_tree(fastq_dir)
    for fi in fastq_files:
        suffix = fi.split(".fastq")[0]
        assert os.path.isfile(suffix + ".bam")
        suffix2 = suffix + "_sorted.bam"
        assert os.path.isfile(suffix2)


# if __name__ == "__main__":
#     test_alginment_workflow_locally()