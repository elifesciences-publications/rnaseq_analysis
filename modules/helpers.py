import configparser
import datetime as dt
import os
import pandas as pd
import re


def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode("utf-8")
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def process_config(config_file="config"):  # tested locally

    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = {}
    for section in config.sections():
        config_dict[section] = {name: value for name, value in config.items(section)}
    return config_dict


def get_second(first):

    pe = re.compile(r'[\W_][0-9][\W_]|[\W_]forward[\W_]')
    matched = pe.search(first).group()
    if "1" in matched:
        new_matched = matched.replace("1", "2")
    else:
        new_matched = matched.replace("forward", "reverse")
    second = first.replace(matched, new_matched)
    return second


def find_files_in_a_tree(folder, file_type='fastq'):

    f_files = []
    for root, dirs, files in os.walk(folder, topdown=False):
        for name in files:
            if name.endswith(file_type):
                f_files.append(os.path.join(root, name))
    return f_files


def set_up_logger():
    # todo set up logger/logging
    print("Set up logger!")


def parse_flagstat(flagstat):

    with open(flagstat, "r") as fh:
        line1 = fh.readline()
        total = int(line1.split()[0])
        fh.readline()
        fh.readline()
        fh.readline()
        line5 = fh.readline()
        mapped = int(line5.split()[0])
        prcnt = mapped/total
        return total, mapped, prcnt


def flagstat_summary(flagstat_dir):
    today = dt.datetime.today().strftime("%Y_%m_%d")
    flagstats = find_files_in_a_tree(flagstat_dir, file_type="flagstat.txt")
    stats = []
    labels = ["Name", "Total", "Mapped", "% Mapped"]
    for fi in flagstats:
        sample_name = os.path.basename(fi).split(".")[0]
        total, mapped, prcnt = parse_flagstat(fi)
        stats.append((sample_name, total, mapped, prcnt))
    df = pd.DataFrame.from_records(stats, index="Name", columns=labels)
    filename = os.path.join(flagstat_dir, todagiy+"_alignment_stats.csv")
    df.to_csv(filename)
    return filename


def process_bedtools_count_output(count_file, count_file_out, param_dict):
    feat = param_dict["feat"]
    count_dict = {}
    with open(count_file, "r") as fh:
        for line in fh:
            if feat not in line:
                continue
            else:
                gene_info = line.split("\t")[8].split(feat+"=")[1].split(";")[0].strip()
                counts = line.split("\t")[9].strip()
                coverage = line.split("\t")[-1].strip()
                count_dict[gene_info] = (counts, coverage)
    with open(count_file_out, "w") as fo:
        for key, val in count_dict.items():
            fo.write("{},{},{}\n".format(key, val[0], val[1]))
    return count_file_out

if __name__ == "__main__":
    flagstat_summary("/Users/annasintsova/git_repos/"
                     "spatial_dynamics_of_gene_expression_in_response_to_T6SS_attack/data/bams")
