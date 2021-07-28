from Bio.Seq import translate, reverse_complement
from Bio import Align
import Bio.Data.CodonTable
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from typing import NamedTuple
import pandas as pd


class Gene(NamedTuple):
    start: int
    end: int
    length: int
    type: str
    trans_table: int
    translation: str
    gene: str
    locus_tag: str
    gc_percentage: float
    strand: int
    positive_percentage: float
    negative_percentage: float
    product: str


def get_list_of_nucleotides(letter):
    nucleotides = ['A', 'C', 'G', 'T']
    nucleotides.remove(letter)
    return nucleotides


def get_count_by_letter(table, amino, nucleotides, codon: str):
    count = 0
    non_stop_codon = 0
    for nucleotide in nucleotides:
        new_codon = codon.format(nucleotide)
        if table[new_codon] == amino:
            count += 1
            non_stop_codon += 1
        elif table[new_codon] != amino and table[new_codon] != '_':
            non_stop_codon += 1
    return count, non_stop_codon


def get_dict_table(table_id: int):
    codons = Bio.Data.CodonTable.generic_by_id[table_id]
    table = codons.forward_table
    stop_codons = codons.stop_codons
    for codon in stop_codons:
        if 'U' not in codon:
            table[codon] = '_'
    for key in list(table):
        if 'U' in key:
            del table[key]
    return table


def count_codon(table, codon):
    synonymous = 0
    non_stop_codon_to_dict = 0
    amino = table[codon]
    letters = [codon[0:1], codon[1:2], codon[2:3]]
    for i, letter in enumerate(letters):
        nucleotides = get_list_of_nucleotides(letter)
        new_letters = letters.copy()
        new_letters[i] = "{}"
        count, non_stop_codon = get_count_by_letter(table, amino, nucleotides, ''.join(new_letters))
        synonymous += count
        non_stop_codon_to_dict += non_stop_codon
    return 3 * synonymous / non_stop_codon_to_dict


def get_synonymous_from_codon(table_id: int):
    amino_acids = {}
    table = get_dict_table(table_id)

    for key in table.keys():
        synonymous = count_codon(table, key)
        amino_acids[key] = synonymous

    print(amino_acids)
    return amino_acids


def shared_genes(df, df2):
    df3 = df[(df.gene != '') & (df.type != 'gene')]
    df4 = df2[(df2.gene != '') & (df2.type != 'gene')]
    first_seq_genes = df3['gene']
    second_seq_genes = df4['gene']
    shared_genes_name = []
    count = 0
    for gene in first_seq_genes:
        if gene in second_seq_genes.values:
            count += 1
            shared_genes_name.append(gene)
    print("Shared genes between the two sequences is: {}".format(count))
    return shared_genes_name


def get_dn_ds_ratio(new_seq, new_seq2, codon_table):
    dn, ds = cal_dn_ds(new_seq, new_seq2, codon_table=codon_table)
    dn_ds = (dn / ds)
    selection = what_is_the_selection(dn_ds)
    return dn_ds, selection


def get_date(df, df2, share_gene):
    df3 = df[(df['gene'] == share_gene) & (df['type'] != 'gene')]
    gene1 = Gene._make(df3.iloc[0, :])
    df4 = df2[(df2['gene'] == share_gene) & (df2['type'] != 'gene')]
    gene2 = Gene._make(df4.iloc[0, :])
    return gene1, gene2


def get_gap_sequences_by_protein_alignment(seq, gene1, seq2, gene2, codon_table):
    gene_seq = get_seq_and_reverse_if_needed(seq, gene1.start, gene1.end, gene1.strand)
    gene_seq2 = get_seq_and_reverse_if_needed(seq2, gene2.start, gene2.end, gene2.strand)

    align_p, second_align_p = get_aligned_proteins(gene_seq, gene_seq2, codon_table)
    if align_p == "error":
        return "error", "error"

    new_seq = insert_gaps(align_p, gene_seq)
    new_seq2 = insert_gaps(second_align_p, gene_seq2)
    new_seq = CodonSeq(str(new_seq))
    new_seq2 = CodonSeq(str(new_seq2))
    return new_seq, new_seq2


def get_aligned_proteins(gene_seq, second_gene_seq, codon_table):
    aligner = Align.PairwiseAligner()
    trans_one = translate(gene_seq, table=codon_table)[:-1]
    trans_two = translate(second_gene_seq, table=codon_table)[:-1]

    try:
        alignments = aligner.align(trans_one, trans_two)
    except:
        return "error", "error"

    splinter = str(alignments[0]).split('\n')
    align_p = splinter[0]
    second_align_p = splinter[2]
    return align_p, second_align_p


def get_seq_and_reverse_if_needed(seq, start, end, strand):
    gene_seq = seq[start:end]
    if strand == -1:
        gene_seq = reverse_complement(gene_seq)
    return gene_seq


def insert_gaps(align_p, gene_seq):
    j = 0
    new_seq = ""
    for char in align_p:
        if char != '-':
            new_seq += gene_seq[j:j + 3]
            j += 3
        else:
            new_seq += '---'
    return new_seq


def what_is_the_selection(dn_ds):
    if dn_ds >= 1.05:
        return "positive."
    elif 1.05 > dn_ds > 0.95:
        return "neutral."
    else:
        return "negative."


def get_sequences(df, df2, shared_genes_names, seq, seq2):
    gene_data = {}
    for i, shared_gene in enumerate(shared_genes_names):

        gene1, gene2 = get_date(df, df2, shared_genes_names[i])
        codon_table = Bio.Data.CodonTable.generic_by_id[int(gene1.trans_table)]

        new_seq, new_seq2 = get_gap_sequences_by_protein_alignment(seq, gene1, seq2, gene2, codon_table)
        if new_seq == "error":
            continue

        dn_ds, selection = get_dn_ds_ratio(new_seq, new_seq2, codon_table)
        gene_data[gene1.gene] = {
            "name": gene1.gene,
            "product": gene1.product,
            "length": gene1.length,
            "translation": gene1.translation,
            "strand": gene1.strand,
            "dn_ds_ratio": dn_ds,
            "The selection is": selection
        }
    df = pd.DataFrame(gene_data, columns=gene_data.keys()).transpose()
    df.to_csv("dn_ds.csv")
