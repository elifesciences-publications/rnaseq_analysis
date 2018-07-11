import sys
sys.path.append('.')
from modules import keggAPI

def test_get_genes():
    gene_file = "/Users/annasintsova/git_repos/code/tests/test_data/2018-06-27-Dienes-line-4-hours-time.csv"

    actual_list = keggAPI.get_genes("pmr", gene_file)

    final_list =   ["pmr:PMI1031",
                    "pmr:PMI0781",
                    "pmr:PMI0348",
                    "pmr:PMI2408",
                    "pmr:PMI1956",
                    "pmr:PMI2662",
                    "pmr:PMI3226",
                    "pmr:PMI3039",
                    "pmr:PMI1804",
                    "pmr:PMI3549",
                    "pmr:PMI1011",
                    "pmr:PMI3550"]
    assert actual_list == final_list