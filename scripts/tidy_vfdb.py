import sys
import json
from os import path

import pandas as pd
from Bio import SeqIO


# blast out columns names
# ["query_name", "subject_name", "identities", "align", "diff_bases", "gap", "query_start",
#  "query_end", "subject_start", "subject_end", "evalue", "bit_score"]

def tidy_blast(blast_out, vfdb_json, predict_nuc, merged_results, outdir):
    blast_df = pd.read_csv(blast_out, sep="\t", header=None,
                           names=[
                               "query_name", "subject_name", "identities", "align", "diff_bases", "gap",
                               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"]
                           )

    vfdb_df = pd.DataFrame.from_dict(json.load(open(vfdb_json, 'r')), orient="index").reset_index(drop=True)
    result_df = blast_df.merge(vfdb_df, left_on="subject_name", right_on="gene_id")
    result_df.to_csv(merged_results, sep="\t", index=False)
    query_names = set(result_df['query_name'])
    records = SeqIO.to_dict(SeqIO.parse(open(predict_nuc, "r"), "fasta"))
    for query_name in query_names:
        with open(path.join(outdir, f"{query_name}.fasta"), "w", encoding='utf-8') as outf:
            outf.write(f">{query_name}\n")
            outf.write(str(records[query_name].seq) + "\n")


if __name__ == '__main__':
    blast_out, vfdb_json, predict_nuc, merged_results, outdir = sys.argv[1:]
    tidy_blast(blast_out, vfdb_json, predict_nuc, merged_results, outdir)
