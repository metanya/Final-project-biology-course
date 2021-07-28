from collections import Counter
from Bio.Seq import translate, reverse_complement
import Bio.Data.CodonTable
import matplotlib.pyplot as plt
import pandas as pd


def count_areas(df) -> {}:
    column = df['type']
    a = Counter(column)
    return dict(a)


def get_areas_mean(df) -> []:
    df2 = df[df['type'] == 'CDS']
    column = df2['length']
    return column.mean()


def print_the_histogram_by_list(column, title):
    a = Counter(column)
    plt.bar(a.keys(), a.values(), color='purple')
    plt.title(title)
    # The categories represent the lengths of the genes and the amount represent how many genes have the same length.
    plt.xlabel('Category')
    plt.ylabel('Amount')
    plt.show()


def print_histograms(df):
    df2 = df[df['type'] == 'CDS']
    print_the_histogram_by_list(df2['length'], 'Lengths of proteins.')
    df2 = df[(df['type'] != 'CDS') & (df['type'] != 'source')]
    print_the_histogram_by_list(df2['length'], 'Lengths of not proteins.')


def get_gc_percentage_mean(df):
    df2 = df[df['type'] == 'CDS']
    column = df2['gc_percentage']
    return column.mean()


def print_gc_percentage_histogram(df):
    df2 = df[df['type'] == 'CDS']
    print_the_histogram_by_list(df2['gc_percentage'], 'Lengths of gc percentage.')


def is_start_codon(seq_part, codon_table, features, start, locus_tag):
    name = str(start) + locus_tag
    if seq_part[0:3] not in codon_table.start_codons and seq_part[1:4] not in codon_table.start_codons \
            and seq_part[2:5] not in codon_table.start_codons:
        features[start] = {"name": name,
                           "cause": "gene not starting with start codon. name in csv is: {} CDS "
                                    "and his locus_tag is: {}. ".format(start, locus_tag)}


def conflict(df, seq):
    df = df[df['type'] == 'CDS']
    start = df['start'].values
    end = df['end'].values
    protein = df['translation'].values
    table = df['trans_table'].values
    strand = df['strand'].values
    locus_tag = df['locus_tag'].values
    features = {}
    for start, end, protein, table, strand, locus_tag in zip(start, end, protein, table, strand, locus_tag):

        codon_table = Bio.Data.CodonTable.generic_by_id[int(table)]
        seq_part = seq[start:end]
        name = str(start) + locus_tag

        if strand == -1:
            seq_part = reverse_complement(seq_part)
        try:
            translation = translate(seq_part, table=codon_table)[:-1]
        except:
            features[start] = {"name": name,
                               "cause": "Positive - couldn't be able to translate. name in csv is: {} CDS and his "
                                        "locus_tag is: {}".format(start, locus_tag)}
            continue
        is_start_codon(seq_part, codon_table, features, start, locus_tag)
        translation = 'M' + translation[1:]
        if translation == protein:
            continue
        elif translation == '' or protein == '':
            features[start] = {"name": name,
                               "cause": "One of the proteins don't contain protein. name in csv is: {} CDS "
                                        "and his locus_tag is: {}. ".format(start, locus_tag)}
        else:
            features[start] = {"name": name,
                               "cause": "Difference in csv protein and translation. name in csv is: {} CDS "
                                        "and his locus_tag is: {}. ".format(start, locus_tag)}

    df = pd.DataFrame(features, columns=features.keys()).transpose()
    df.to_csv("gene_exceptions.csv")


def get_positive_negative_mean(df) -> []:
    df2 = df[df['type'] == 'CDS']
    df2 = df2[df2['amino_acid_positive_percentage'] != -1]
    df2 = df2[df2['amino_acid_negative_percentage'] != -1]
    positive = df2['amino_acid_positive_percentage']
    print_the_histogram_by_list(positive, 'Lengths of positive amino acids percentage.')
    negative = df2['amino_acid_negative_percentage']
    print_the_histogram_by_list(negative, 'Lengths of negative amino acids percentage.')
    return positive.mean(), negative.mean()


class PartA:

    @staticmethod
    def first_question(df):
        areas = count_areas(df)
        print("Dictionary of the genome areas: {}.".format(areas))

    @staticmethod
    def second_question(df):
        mean = get_areas_mean(df)
        print("Average of proteins length is: {}.".format(mean))
        print_histograms(df)

    @staticmethod
    def third_question(df):
        mean = get_gc_percentage_mean(df)
        print("Average gc percentage in the sequence is: {} and average gc percentage of proteins is: {}".format(
            df.loc['0 source', 'gc_percentage'], mean))
        print("The outcome is weird to me because we said gc percentage help the genes and make them more stable.")
        print("I thought the proteins should contain great amount of G and C cause they really important.")
        print("Therefore I think The small GC content show on some capabilities of the bacteria.")

        print_gc_percentage_histogram(df)

        df3 = df[(df['type'] != 'gene')]
        print("5 richest genes with gc percentage is:\n {}.".format(
            df3[['gc_percentage', 'gene']].sort_values(by=['gc_percentage'], ascending=False).head(5)))
        print("5 poorest genes with gc percentage is:\n {}.".format(
            df3[['gc_percentage', 'gene']].sort_values(by=['gc_percentage'], ascending=True).head(5)))

    @staticmethod
    def fourth_question(df, seq):
        conflict(df, seq)
        positive, negative = get_positive_negative_mean(df)
        print("The mean of the positive amino acids in the proteins is: {} and the negative is: {}.".format(positive,
                                                                                                            negative))
