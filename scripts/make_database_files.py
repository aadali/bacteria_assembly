import os
import re
import json
import sys

from os import path
from collections import defaultdict

from Bio import SeqIO

"""
THIS SCRIPT WILL NOT BE USED IN ASSEMBLY PIPELIE

There some functions were used to make some pubmlst and vfdb database files. 
Such as scheme.txt, locus_vs_scheme.json, st_alleles.json of pubmlst
Such as vfdb.json of vfdb
"""


def get_vfdb_info(setA_nt, setA_pro):
    """
    get vfdb.json to store all the info of virulence genes or protein
    :param setA_nt: the core virulence genes from VFDB
    :param setA_pro: the core virulence protein from VFDB
    :return:
    """
    nuc_records = SeqIO.parse(open(setA_nt, "r"), "fasta")
    pro_records = SeqIO.parse(open(setA_pro, "r"), "fasta")
    seq_dict = {}
    for nuc_record in nuc_records:
        seq_dict[nuc_record.id] = [str(nuc_record.seq)]

    for pro_record in pro_records:
        seq_dict[pro_record.id].append(str(pro_record.seq))

    pat = re.compile("(VFG\d+\(gb.+?\)) \((.+)\) .+ \[(.+) \((VF\d+)\) - (.+) \((VFC\d+)\)\] \[(.+)\]", re.IGNORECASE)
    records = SeqIO.parse(open(setA_pro, 'r', encoding='utf-8'), "fasta")
    vfg = {}
    for record in records:
        mat = pat.search(record.description)
        try:
            gene_id = mat.group(1)
            gene_name = mat.group(2)
            vf_name = mat.group(3)
            vf_id = mat.group(4)
            vfc_name = mat.group(5)
            vfc_id = mat.group(6)
            bacteria = mat.group(7)
            vfg[gene_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "vf_name": vf_name,
                "vf_id": vf_id,
                "vfc_name": vfc_name,
                "vfc_id": vfc_id,
                "bacteria": bacteria,
                "nt_seq": seq_dict[gene_id][0],
                'pro_seq': seq_dict[gene_id][1]
            }
        except (IndexError, AttributeError) as e:
            print(e, file=sys.stderr)
            print(f"{record.name} error when construct VFG dict", file=sys.stderr)
            continue
    json.dump(vfg, open("vfdb.json", 'w'))


def download(xmlfile):
    """
    parse the dbases.xml from https://pubmlst.org/data , generate the download.sh to download the profile.tsv and
    each alleles fasta file
    After the download.sh generated, do `bash download.sh`
    """
    # xmlfile = "../database/pubmlst/mlst_dbases.xml"
    f = open(xmlfile, "r")
    xml = f.read()
    mat = re.findall("<species>(.+?)</mlst>", xml, flags=re.DOTALL)
    shell = open("download.sh", "w", encoding='utf-8')
    for record in mat:
        specie = re.findall("(.+)\n<mlst>", record)
        specie = specie[0].replace(" ", "_").replace("(", "_").replace(")", "_").replace("/", "#")
        profies_url = re.findall("<profiles>\n<url>(.+?)</url>", record)
        locus = re.findall("<locus>(.+?)\n<url>", record)
        locus_url = re.findall("<locus>.+?\n<url>(.+?)</url>", record)
        shell.write(f"mkdir {specie}\n")
        shell.write(f"wget -O ./{specie}/profile.tsv {profies_url[0]}\n")
        for each_locus, each_locus_url in zip(locus, locus_url):
            shell.write(f"wget -O ./{specie}/{each_locus}.fasta {each_locus_url}\n")
            shell.write("\n\n")
        assert len(locus) == len(locus_url)


def get_scheme(xmlfile):
    """
    get tsv file: column1 is the scheme name, column2 is the locuses belong to the scheme
    :param xmlfile: dbases.xml from https://pubmlst.org/data
    :return:
    """
    # xml_file = "/home/a/big/ycq/projects/assembly/database/pubmlst/mlst_dbases.xml"

    file = open(xmlfile, 'r', encoding='utf-8')
    xml = file.read()

    mat = re.findall("<species>(.+?)</mlst>", xml, flags=re.DOTALL)
    schemes = {}
    for record in mat:
        specie = re.findall("(.+)\n<mlst>", record)
        specie = specie[0].replace(" ", "_").replace("(", "_").replace(")", "_").replace("/", "#")
        profies_url = re.findall("<profiles>\n<url>(.+?)</url>", record)
        locus = re.findall("<locus>(.+?)\n<url>", record)
        locus_url = re.findall("<locus>.+?\n<url>(.+?)</url>", record)
        schemes[specie] = locus

    with open("schemes.txt", 'w', encoding='utf-8') as outfile:
        for key, value in schemes.items():
            outfile.write(key + "\t")
            outfile.write(",".join(value) + '\n')


def get_st_alleles(database_dir, scheme_file):
    """
    get the json file. All the info about alleles vs st of each scheme were collected.
    Each scheme record in the json contain info:
        1) "locuses": a list containing all locuses of this scheme
        2) "alleles": key: all the different locus alleles group
                      value: st and something
    :param database_dir: the path where download.sh worked. Only the scheme directory allowed existing in this directory
    :param scheme_file: the outfile of get_scheme(xmlfile)
    :return:
    """
    dirs = os.listdir(database_dir)
    dirs = list(filter(lambda x: path.isdir(path.join(database_dir, x)), dirs))
    bases_names = [dir for dir in dirs]
    d = {}
    with open(scheme_file, 'r') as infile:
        for line in infile:
            new_line = line.strip()
            key, values = new_line.split("\t")
            d[key] = {}
            d[key]['locuses'] = values.split(",")
            d[key]["alleles"] = {}

    for dir in bases_names:
        locus_numbers = len(d[dir]['locuses'])
        with open(path.join(f"{database_dir}/{dir}", "profile.tsv"), 'r') as infile:
            header = next(infile).strip("\n").split("\t")
            if "clonal_complex" in header:
                clonal_complex_idx = header.index('clonal_complex')
            else:
                clonal_complex_idx = None

            if "species" in header:
                species_id = header.index("species")
            else:
                species_id = None

            locus = header[1:1 + locus_numbers]
            # print(locus, len(locus))
            for line in infile:
                new_line = line.strip("\n")
                st_alleles = new_line.split("\t")

                st = st_alleles[0]
                alleles = st_alleles[1:1 + locus_numbers]
                locus_allele = [f"{loci}@{allele}" for loci, allele in zip(locus, alleles)]

                st_key = "|".join(list(sorted(locus_allele)))
                d[dir]["alleles"][st_key] = {}
                d[dir]["alleles"][st_key]['st'] = st
                if not clonal_complex_idx is None:
                    d[dir]["alleles"][st_key]['clonal_complex'] = st_alleles[clonal_complex_idx]
                if not species_id is None:
                    d[dir]["alleles"][st_key]['clonal_complex'] = st_alleles[species_id]
    json.dump(d, open("st_alleles.json", 'w'))


def get_locus_vs_scheme(scheme_file):
    """
    get the json file: each locus belongs to which schemes
    :param scheme_file: the outfile of get_scheme(xmlfile)
    :return:
    """
    d = defaultdict(list)
    with open(scheme_file, 'r') as infile:
        for line in infile:
            scheme, locuses = line.strip().split("\t")
            locuses = locuses.split(",")
            for locus in locuses:
                if scheme not in d[locus]:
                    d[locus].append(scheme)

    json.dump(d, open("locus_vs_scheme.json", 'w'))
