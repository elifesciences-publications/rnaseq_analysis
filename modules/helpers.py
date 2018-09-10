import configparser
import os
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