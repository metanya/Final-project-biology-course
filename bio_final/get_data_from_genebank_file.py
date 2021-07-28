from Bio import SeqIO
import pandas as pd


def get_protein_gc_percentage(feature, seq):
    seq = seq[feature.location.start.position: feature.location.end.position]
    counter_g = seq.count('G')
    counter_c = seq.count('C')
    gc_percentage = ((counter_c + counter_g) / len(seq)) * 100
    return gc_percentage


def positive_negative_amino_acids(translation):
    trans_length = len(translation)
    counter_e = translation.count('E')
    counter_d = translation.count('D')
    counter_k = translation.count('K')
    counter_h = translation.count('H')
    counter_r = translation.count('R')
    positive = ((counter_r + counter_h + counter_k) / trans_length) * 100
    negative = ((counter_e + counter_d) / trans_length) * 100
    return positive, negative


class GetDataFromGenebankFile:

    @staticmethod
    def get_data_frame(gene_bank_file: str, csv_name: str):
        features = {}
        with open(gene_bank_file, "r") as input_handle:
            for record_gb in SeqIO.parse(input_handle, "genbank"):
                seq = record_gb.seq.upper()
                for feature in record_gb.features:
                    trans_table = ""
                    translation = ""
                    gene = ""
                    locus_tag = ""
                    positive = -1
                    negative = -1
                    product = ""
                    if 'transl_table' in feature.qualifiers:
                        trans_table = feature.qualifiers['transl_table'][0]
                    if 'translation' in feature.qualifiers:
                        translation = feature.qualifiers['translation'][0]
                        positive, negative = positive_negative_amino_acids(translation)
                    if 'gene' in feature.qualifiers:
                        gene = feature.qualifiers['gene'][0]
                    if 'locus_tag' in feature.qualifiers:
                        locus_tag = feature.qualifiers['locus_tag'][0]
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]

                    gc_percentage = get_protein_gc_percentage(feature, seq)
                    name = str(feature.location.start) + " " + feature.type
                    features[name] = {"start": feature.location.start,
                                      "end": feature.location.end,
                                      "length": feature.location.end - feature.location.start,
                                      "type": feature.type,
                                      "trans_table": trans_table,
                                      "translation": translation,
                                      "gene": gene,
                                      "locus_tag": locus_tag,
                                      "gc_percentage": gc_percentage,
                                      "strand": feature.location.strand,
                                      "amino_acid_positive_percentage": positive,
                                      "amino_acid_negative_percentage": negative,
                                      "product": product}
        df = pd.DataFrame(features, columns=features.keys()).transpose()
        df.to_csv(csv_name)
        return df, seq
