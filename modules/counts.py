

def count_reads(file, out_dir, ref, config,
                form="bam", order="pos", attr = "ID",
                mode="union", stranded="yes",
                feature = "CDS"):
    """


    :param files:
    :param out_dir:
    :param ref: what is a ref? gff file?
    :param config:
    :param form:
    :param order:
    :param mode:
    :param stranded:
    :param feature:
    :return:
    """

    # htseq_path = config.get("HTSEQ", "bin")
    scripts = []
    if os.path.isdir(ref):
        annot_files = [os.path.join(ref, r) for r in os.listdir(ref)]
    else:
        annot_files = [ref]
    for fh in files:
        prefix = os.path.basename(fh).split("_")[0]
        if len(annot_files) > 1:
            annot = [an for an in annot_files if prefix in an][0]
        else:
            annot = annot_files[0]
        count_file = os.path.join(out_dir, (os.path.basename(fh).split(".")[0] + "_counts"))
        script = "htseq-count -f {0} -r {1} -m {2}" \
                  " -s {3} -t {4} -i {8} {5} {6} > {7}\n".format(form, order,
                                                              mode, stranded, feature, fh,
                                                              annot, count_file, attr)
        scripts.append((os.path.basename(count_file), script))

    return scripts
