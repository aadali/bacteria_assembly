import sys
import json

import pandas as pd


# the below lines is the headers of blast out and what does it mean
# sseqid slen length nident pident qseqid qstart qend qseq sstrand
#  sseqid means Subject Seq-id
#    slen means Subject sequence length
#  length means Alignment length
#  nident means Number of identical matches
#  pident means Percentage of identical matches
#  qseqid means Query Seq-id
#  qstart means Start of alignment in query
#   qend  means End of alignment in query
#   qseq  means means Aligned part of query sequence
# sstrand means Subject Strand

class AlignAllele(object):
    """
    Struct to store a locus align info
    """

    def __init__(self, locus, allele, qseqid, qstart, qend, pident):
        self.locus = locus
        self.allele = allele
        self.qseqid = qseqid
        self.qstart = qstart
        self.qend = qend
        self.pident = pident

    def __repr__(self):
        info = dict(
            locus=self.locus,
            allele=self.allele,
            qseqid=self.qseqid,
            qstart=self.qstart,
            qend=self.qend,
            pident=self.pident
        )
        return f"{info}"


class Parser(object):
    """
    input pubmlst blast out, output the mlst result
    """

    def __init__(self, blast_out, locus_scheme_js, st_alleles_js):
        """
        init
        :param blast_out: the blast out file by blastn
        :param locus_scheme_js: the locus_vs_scheme.json which was prepared in advance
        :param st_alleles_js: the st_alleles.json which was prepared in advance
        """
        self.blast_df = pd.read_csv(
            blast_out, sep="\t", header=None,
            names=["sseqid", "slen", "length", "nident", "pident", "qseqid", "qstart", "qend", "qseq", "sstrand"],
            dtype={"sseqid": str,
                   "slen": int,
                   "length": int,
                   "nident": int,
                   "pindent": float,
                   "qseqid": str,
                   'qstart': int,
                   "qend": int,
                   "qseq": str,
                   "sstrand": str}
        )
        self.locus_scheme = json.load(open(locus_scheme_js, 'r'))
        self.st_alleles = json.load(open(st_alleles_js, 'r'))

    def _get_blast_schemes(self):
        """
        goouped_by locus, and get the best align of each locus, drop the align whose pident less than 95
        :return: a sub dataframe
        """
        locus_allele = self.blast_df['sseqid'].str.rsplit("_", 1, expand=True)
        locus_allele.columns = ['locus', 'allele']
        self.blast_df.insert(loc=0, column="locus", value=locus_allele['locus'])
        self.blast_df.insert(loc=1, column="allele", value=locus_allele['allele'])
        grouped = self.blast_df.groupby("locus")
        good_blast_records = grouped.head(1).sort_values(by="pident", ascending=False).query("pident > 95")
        return good_blast_records

    def _init_scores(self):
        """
        construct a directory to store the scheme info from blast out
        :return: the candidate schemes
        """
        good_blast_records = self._get_blast_schemes()
        schemes = [y for x in good_blast_records['locus'] for y in self.locus_scheme[x]]
        scores = {}
        for scheme in set(schemes):
            scores[scheme] = {}  # construcn a empty dict
            scores[scheme]['total_score'] = 0  # each locus' pident of this scheme will be added
            scores[scheme]['alleles'] = []  # to store real locus thaf found from blast out
            count_locus = len(self.st_alleles[scheme]['locuses'])
            scores[scheme]['locus_number'] = count_locus  # how many locus this scheme has
            scores[scheme]['max_score'] = 100 * count_locus  # if all the locus is exact match, the score will be 100
        return good_blast_records, scores

    def write_mlst(self, bad_align_num=1):
        good_blast_records, scores = self._init_scores()
        for row in good_blast_records.itertuples(index=False):
            locus = row.locus
            allele = row.allele
            pident = row.pident
            qseqid = row.qseqid
            qstart = row.qstart
            qend = row.qend
            for scheme in self.locus_scheme[locus]:
                scores[scheme]['total_score'] = scores[scheme]['total_score'] + pident
                scores[scheme]['alleles'].append(
                    AlignAllele(locus=locus, allele=allele, pident=pident, qseqid=qseqid, qstart=qstart, qend=qend)
                )
                scores[scheme]['score'] = scores[scheme]['total_score'] / scores[scheme]['max_score'] * 100

        # sorted by the score,
        sorted_scores = list(sorted(scores.items(), key=lambda x: x[1]['score'], reverse=True))
        first = sorted_scores[0]
        st_find = False
        for idx, (scheme, scheme_aligns) in enumerate(sorted_scores):
            score = scheme_aligns['score']
            locuses = [x.locus for x in scheme_aligns['alleles']]
            alleles = [x.allele for x in scheme_aligns['alleles']]
            pidents = [x.pident for x in scheme_aligns['alleles']]
            sorted_locuses = list(
                sorted(locuses))  # sorted the locus to make the key of the alleles of the scheme from st_alleles.json
            sorted_alleles = [alleles[locuses.index(each)] for each in
                              sorted_locuses]  # sorted the allele to make the key of the alleles of the scheme from st_alleles.json
            # sroted_pidents = [pidents[locuses.index(each)] for each in sorted_locuses]
            expect_locuses_numbers = len(self.st_alleles[scheme]['locuses'])
            real_bad_pidents = list(filter(lambda x: x < 99, pidents))

            # if too many bad locus align was found, this scheme will be dropped
            if (real_bad_align_num := len(real_bad_pidents)) > bad_align_num:
                print(f"{scheme} was filted because too many bad locus aligns: {real_bad_align_num} found",
                      file=sys.stderr)
                continue

            # if too few locsus was found, this scheme will be dropped
            if (real_locuses_num := len(sorted_locuses)) < expect_locuses_numbers:
                print(
                    f"{scheme} was filted because not enough locuses, {real_locuses_num} found, but expected {expect_locuses_numbers}",
                    file=sys.stderr)
                continue

            # all the locus align of first scheme is exact, it's a good result
            # Maybe, there will be many scheme found, such as Escherichia_coli# and Escherichia_coli#2
            if first[1]['score'] == 100 and score == 100:
                st_find=True
                allele_key = "|".join([f"{x}@{y}" for x, y in zip(sorted_locuses, sorted_alleles)])
                st = self.st_alleles[scheme]["alleles"][allele_key]['st']
                print("\t".join(["Scheme", "Locus", "Allele", "Identity", "Contig", "Start", "End", "ST"]))
                for align_allele in scheme_aligns["alleles"]:
                    print("\t".join(
                        [scheme, align_allele.locus, align_allele.allele, str(align_allele.pident), align_allele.qseqid,
                         str(align_allele.qstart),
                         str(align_allele.qend), str(st)]))
            else:
                # for the midding good and midding bad align, there will be not ST. "-" was used to instead of st number
                # And the specified locus will be marked by the prefix: "~"
                st_find =True
                print("\t".join(["Scheme", "Locus", "Allele", "Identity", "Contig", "Start", "End", "ST"]))
                for align_allele in scheme_aligns["alleles"]:
                    print("\t".join([
                        scheme,
                        align_allele.locus,
                        align_allele.allele if align_allele.pident == 100 else f"~{align_allele.allele}",
                        str(align_allele.pident),
                        align_allele.qseqid,
                        str(align_allele.qstart),
                        str(align_allele.qend),
                        "-"
                    ]))
        if not st_find:
            print("###No ST found")


if __name__ == '__main__':
    # js_file1 = "/home/a/big/ycq/projects/assembly/database/pubmlst/locus_vs_scheme.json"
    # js_file2 = "/home/a/big/ycq/projects/assembly/database/pubmlst/st_alleles.json"
    # blast_out = "/home/a/Desktop/blast_out.tsv"
    # blast_out = "/home/a/Desktop/anyang-bar01.mlst.tsv"
    blast_out, locus_scheme, st_alleles = sys.argv[1:]
    p = Parser(blast_out=blast_out, locus_scheme_js=locus_scheme, st_alleles_js=st_alleles)
    p.write_mlst()
