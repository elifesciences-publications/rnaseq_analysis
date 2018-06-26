def match_locus_tags(gff, matched_tags_file):
    with open(matched_tags_file, "w") as fo:
        with open(gff, "r") as fi:
            for line in fi:
                if not "locus_tag" in line:
                    continue
                else:
                    #print(line)
                    if "old_locus_tag" in line:
                        old_tag = line.split("old_locus_tag=")[1].split(";")[0].strip()
                    else:
                        old_tag = "NA"
                    new_tag = line.split("old_locus_tag=")[0].split("locus_tag=")[1].split(";")[0].strip()
                    fo.write("{},{}\n".format(old_tag, new_tag))



if __name__ == "__main__":
    gff = "/Users/annasintsova/git_repos/proteus/data/annotations/GCF_000069965.1_ASM6996v1_genomic.gff"
    matched_tags_file = "/Users/annasintsova/git_repos/proteus/data/annotations/new_locus_tags.csv"
    match_locus_tags(gff, matched_tags_file)