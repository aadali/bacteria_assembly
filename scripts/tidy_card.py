import sys
from os import path
import pandas as pd


def tidy_diamond(diamond_out, tidy_file, aro_index_file, outdir):
    aro_df = pd.read_csv(aro_index_file, sep="\t", names=
    ['ARO_Accession', 'CVTERM_ID', 'ModelSequenceID', 'ModelID', 'ModelName', 'ARO_Name',
     'ProteinAccession', 'DNA_Accession', 'AMR_GeneFamily', 'DrugClass', 'ResistanceMethanism',
     'CARD_ShortName'])
    with open(diamond_out, 'r', encoding='utf-8') as infile, open(tidy_file, 'w', encoding='utf-8') as outfile:
        for line in infile:
            if line.startswith("ORF_ID"):
                continue
            (
            orf_id, contig, start, stop, strand, cut_off, _, _, aro_name, identities, aro, model_type, _, _, drug_class,
            resistance_methanism, amr_gene_family, dna_seq, _, _, percentage_len_of_ref_seq, _, _, _,
            _) = line.strip("\n").split("\t")
            contig = contig.strip()
            aro_accession = f"ARO:{aro}"
            sub_aro_df = aro_df.query("ARO_Accession == @aro_accession").reset_index(drop=True)
            if sub_aro_df.empty:
                print(f"Couldn't find record for {aro}", file=sys.stderr)
                exit(2)
            if sub_aro_df.shape[0] > 1:
                print(f"Multi records found for {aro}", file=sys.stderr)
                exit(2)
            ref_dna_accession = sub_aro_df['DNA_Accession'][0]
            ref_protein_accession = sub_aro_df['ProteinAccession'][0]
            with open(path.join(outdir, f"{contig}.fasta"), 'w') as contig_file:
                contig_file.write(f">{contig}\n")
                contig_file.write(dna_seq + '\n')
            new_line = [orf_id, contig, start, stop, strand, cut_off, aro_name, identities, aro, model_type, drug_class,
                        resistance_methanism, amr_gene_family, percentage_len_of_ref_seq, ref_dna_accession,
                        ref_protein_accession]
            header = ["orf_id", 'contig', 'start', 'stop', 'strand', 'cut_off', 'identities', 'aro', 'model_type',
                      'drug_class', 'resistance_methanism', 'amr_gene_family', 'percentage_len_of_ref_seq',
                      'ref_dna_accession', 'ref_protein_accession']
            outfile.write("#orf_id: Open Reading Frame identifier (internal to RGI)\n")
            outfile.write("#contig: Source Sequence\n")
            outfile.write("#start: Start co-ordinate of ORF\n")
            outfile.write("#stop: End co-ordinate of ORF\n")
            outfile.write("#strand: Strand of ORF\n")
            outfile.write("#cut_off: RGI Detection Paradigm (Perfect, Strict, Loose)\n")
            outfile.write("#aro_name: ARO term of top hit in CARD\n")
            outfile.write("#identities: Percent identity of match to top hit in CARD\n")
            outfile.write("#aro: ARO accession of match to top hit in CARD\n")
            outfile.write("#model_type: CARD detection model type\n")
            outfile.write("#drug_class: ARO Categorization\n")
            outfile.write("#resistance_methanism: ARO Categorization\n")
            outfile.write("#amr_gene_family: ARO Categorization\n")
            outfile.write("#percentage_len_of_ref_seq: (length of ORF protein / length of CARD reference protein)\n")
            outfile.write("#ref_dna_accession: reference dna accession number of CARD")
            outfile.write("#ref_protein_accession: protein dna accession number of CARD")
            outfile.write("\t".join(header) + '\n')
            outfile.write("\t".join(new_line) + '\n')


if __name__ == '__main__':
    diamond_out, tidy_file, aro_index_file, outdir = sys.argv[1:]
    tidy_diamond(diamond_out, tidy_file, aro_index_file, outdir)
