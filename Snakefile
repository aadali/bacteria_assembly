import os
import re
import gzip
import sys
import yaml
from os import path
from subprocess import run

from scripts.softwares import softwares

SNAKEDIR = path.dirname(path.abspath(workflow.snakefile))
CARD_DB = path.join(SNAKEDIR,"database/card/localDB")
ARO_INDEX = path.join(SNAKEDIR,"database/card/aro_index.tsv")
VFDB = path.join(SNAKEDIR,"database/vfdb/vfdb")
VFDB_JSON = path.join(SNAKEDIR,"database/vfdb/vfdb.json")
PUBMLIST_DB = path.join(SNAKEDIR,"database/pubmlst_alleles")
LOCSU_VS_SCHEME = path.join(SNAKEDIR,"database/pubmlst/locus_vs_scheme.json")
ST_ALLEES = path.join(SNAKEDIR,"database/pubmlst/st_alleles.json")
configfile: path.join(SNAKEDIR,"config.yaml")

# init the other parms start...
INPUT_LONG = config.get("input_long")
INPUT_SHORT = config.get("input_short")
MIN_LEN = int(config.get("min_read_len",1000))
MIN_READ_QUAL = int(config.get("min_read_qual",10))
MODEL = config.get("model","r941_min_hac_g507")
ANALYSIS_NAME = config.get("analysis_name","test001")
DATA_TYPE = config.geet("data_type","ont")
TOOL = config.get("tool","flye")
GENOME_SIZE = config.get("genome_size","")
GENOME_SIZE = f"--genome-size {GENOME_SIZE}" if GENOME_SIZE else ""

THREADS = 8


rule ont_no_ref:
    input:
        a=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}.qc.summary.txt",
        b=f"{ANALYSIS_NAME}/assembly/{ANALYSIS_NAME}.consensus.fasta",
        tidy_card_out=f"{ANALYSIS_NAME}/card/{ANALYSIS_NAME}.card.tsv",
        tidy_vfdb_out=f"{ANALYSIS_NAME}/vfdb/{ANALYSIS_NAME}.vfdb.tsv"

rule ont_ref:
    input:
        a=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}.qc.summary.txt",
        b=f"{ANALYSIS_NAME}/assembly/{ANALYSIS_NAME}.consensus.fasta",
        c=f"{ANALYSIS_NAME}/variants/{ANALYSIS_NAME}.vcf",
        d=f"{ANALYSIS_NAME}/card/{ANALYSIS_NAME}.card.tsv",
        e=f"{ANALYSIS_NAME}/vfdb/{ANALYSIS_NAME}.vfdb.tsv"

rule hybrida:
    input:
        b=f"{ANALYSIS_NAME}/assembly/{ANALYSIS_NAME}.consensus.fasta",


rule fastp:
    input:
        raw_fastq=INPUT_LONG
    output:
        clean_fastq=f"{ANALYSIS_NAME}/clean_data/{ANALYSIS_NAME}.clean.fastq.gz",
        js=f"{ANALYSIS_NAME}/clean_data/{ANALYSIS_NAME}.json",
        html=f"{ANALYSIS_NAME}/clean_data/{ANALYSIS_NAME}.html"
    params:
        length_required=MIN_LEN,
        average_qual=MIN_READ_QUAL
    shell:
        f"{softwares['fastp']} -i {{input.raw_fastq}} -o {{output.clean_fastq}} --length_required {{params.length_required}} "
        f"--length_limit 0 -q 0 -w 4 --average_qual {{params.average_qual}} --html {{output.html}} "
        f"--json {{output.js}}"

rule nanoplot:
    input:
        raw_fastq=INPUT_LONG,
        clean_fastq=rules.fastp.output.clean_fastq
    output:
        outdir=directory(f"{ANALYSIS_NAME}/nanoplot"),
        qc_summary=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}.qc.summary.txt",
        raw_fig=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}_raw.LengthvsQualityScatterPlot_dot.png",
        clean_fig=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}_clean.LengthvsQualityScatterPlot_dot.png",
        raw_stats_file=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}_raw.NanoStats.txt",
        clean_stats_file=f"{ANALYSIS_NAME}/nanoplot/{ANALYSIS_NAME}_clean.NanoStats.txt",
        finish=temporary(f".{ANALYSIS_NAME}.nanoplot.finish")
    shell:
        f"{softwares['python']} {SNAKEDIR}/scripts/nanoplot.py {{input.raw_fastq}} {{output.raw_fig}} {{output.raw_stats_file}};\n "
        f"{softwares['python']} {SNAKEDIR}/scripts/nanoplot.py {{input.clean_fastq}} {{output.clean_fig}} {{output.clean_stats_file}};\n "
        f"paste {{output.raw_stats_file}} {{output.clean_stats_file}} | cut -f 1,2,4 | sed '1cMetrics\\traw\\tclean' > {{output.qc_summary}};\n "
        f"touch {{output.finish}}"

########################### For ont reference
rule assembly:
    input:
        fastq=rules.fastp.output.clean_fastq
    output:
        assembly_out=directiry(f"{ANALYSIS_NAME}/assembly"),
        draft_fasta=f"{ANALYSIS_NAME}/assembly/assembly_draft.fasta"
    threads:
        THREAS
    params:
        genome_size=GENOME_SIZE
    shell:
        f"{softwares['flye']} --nano-raw {{input.fastq}} -o {{output.assembly_out}} -t {{threads}} --min-overlap 2000 -i 2 {{params.genome_size}};\n"
        f"mv {{output.assembly_out}}/assembly.fasta {{output.draft_fasta}};\n"

rule medaka_polish:
    input:
        draft_fasta=rules.assembly.output.draft_fasta,
        fastq=rules.fastp.output.clean_fastq
    output:
        polish_fasta=f"{ANALYSIS_NAME}/assembly/consensus.fasta"
    threads:
        4
    shell:
        f"{softwares['medaka_consensus']} -i {{input.fastq}} -d {{input.draft_fasta}} -o {ANALYSIS_NAME}/assembly -f -t  {{threads}}"

rule predict_genes:
    input:
        polish_fasta=rules.medaka_polish.output.polish_fasta
    output:
        predict_protein=f"{ANALYSIS_NAME}/prodigal/{ANALYSIS_NAME}.prediction.proteion.fasta",
        predict_nuc=f"{ANALYSIS_NAME}/prodigal/{ANALYSIS_NAME}.prediction.nuc.fasta"
    shell:
        f"{softwares['prodigal']}  -a {{output.predict_protein}} -d {{output.predict_nuc}} -i {{input.polish_fasta}} -d /dev/null"

rule card_anno:
    input:
        card_database=CARD_DB,
        predict_nuc=rules.predict_genes.output.predict_nuc
    output:
        card_out=f"{ANALYSIS_NAME}/card/{ANALYSIS_NAME}.txt",
    threads:
        2
    shell:
        f"ln -s {{input.card_database}} localDB;\n"
        f"{softwares['rgi']} main -i {{input.fasta}} --local --exclude_nudge -a DIAMOND -t contig -o {ANALYSIS_NAME}/card/{ANALYSIS_NAME} --clean -n {{threads}}"

rule tidy_card_results:
    input:
        card_out=rules.card_anno.output.card_out,
        aro_index=ARO_INDEX,
        outdir=directory(f"{ANALYSIS_NAME}/card/resistance_genes")
    output:
        tidy_card_out=f"{ANALYSIS_NAME}/card/{ANALYSIS_NAME}.card.tsv"
    shell:
        f"{softwares['python']} {SNAKEDIR}/scripts/tidy_card.py {{input.card_out}}"


rule vfdb_anno:
    input:
        vfdb=VFDB,
        predict_nuc=rules.predict_genes.output.predict_nuc
    output:
        vfdb_out=f"{ANALYSIS_NAME}/vfdb/{ANALYSIS_NAME}.vfdb_blast.out"
    shell:
        f"{softwares['blastn']} -query {{input.predict_nuc}} -db {{input.vfdb}} -outfmt 6 -out {{output.vfdb_out}}"

rule tidy_vfdb_results:
    input:
        blast_out=rules.vfdb_anno.output.vfdb_out,
        vfdb_info=VFDB_JSON,
        predict_nuc=rules.predict_genes.output.predict_nuc,
        outdir=directory(f"{ANALYSIS_NAME}/vfdb/virulence_genes")
    output:
        tidy_vfdb_out=f"{ANALYSIS_NAME}/vfdb/{ANALYSIS_NAME}.vfdb.tsv"
    shell:
        f"{softwares['python']} {SNAKEDIR}/scripts/tydy_vfdb.py {{input.blast_out}} {{input.vfdb_info}} {{input.predict_nuc}} {{output.tidy_vfdb_out}} {{input.outdir}}"


rule blast_pubmlst:
    input:
        pubmlst_database=PUBMLIST_DB,
        query=rules.medaka_polish.output.polish_fasta
    output:
        blast_out=f"{ANALYSIS_NAME}/pubmlst/{ANALYSIS_NAME}.tsv"
    shell:
        f"{softwares['blastn']} -query {{input.query}} -db {{input.pubmlst_database}} -max_target_seqs 100000 -outfmt \"6 sseqid slen length nident qseqid qstart qend qseq sstrand \"> {{output.blast_out}}"

rule get_mlst:
    input:
        blast_out=rules.blast_pubmlst.output.blast_out,
        locus_vs_scheme=LOCSU_VS_SCHEME,
        st_alleles=ST_ALLEES
    output:
        mlst_result=f"{ANALYSIS_NAME}/pubmlst/{ANALYSIS_NAME}.mlst.result.tsv"
    shell:
        f"{softwares['python']} {SNAKEDIR}/scripts/get_st_from_blast.py {{input.blast_out}} {{input.locus_vs_scheme}} {{input.st_alleles}} > {{output.mlst_result}}"


################## For ont with reference
rule map:
    input:
        clean_data=rules.fastp.output.clean_fastq