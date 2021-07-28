import os
from part_a import PartA
from part_b import shared_genes, get_sequences, get_synonymous_from_codon
from get_data_from_genebank_file import GetDataFromGenebankFile


def get_file_and_assert(file_path: str) -> str:
    gene_bank_file = file_path
    assert (os.path.exists(gene_bank_file))
    return gene_bank_file


def main():
    path = get_file_and_assert('BS3610.gb')
    path2 = get_file_and_assert('AP006627.gb')
    df, seq = GetDataFromGenebankFile.get_data_frame(path, "part_a.csv")
    df2, seq2 = GetDataFromGenebankFile.get_data_frame(path2, "part_b.csv")

    PartA.first_question(df)
    PartA.second_question(df)
    PartA.third_question(df)
    PartA.fourth_question(df, seq)

    # Second part
    get_synonymous_from_codon(11)
    shared_genes_names = shared_genes(df, df2)
    get_sequences(df, df2, shared_genes_names, seq, seq2)


if __name__ == "__main__":
    main()
