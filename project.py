"""
author: [Tamim Hussein]
studentnumber: [2502572]
"""

import sys
from pathlib import Path
from ingest_data import read_csv, clean_data, generate_sequences
from graph_builder import construct_graph, is_valid_graph, plot_graph
from dna_constructor import construct_dna_sequence
from dna_constructor import construct_all_dna, save_all_dna
from dna_constructor import align_segments


def save_output(s: str, filename: str) -> None:
    with open(filename, "w") as f:
        f.write(s)


def main():

    input_file = sys.argv[1]

    input_path = Path(input_file)

    dna_id, k, df = read_csv(input_file)

    graph_filename = input_path.stem + ".png"
    output_filename = input_path.stem + ".txt"
    output_all_filename = f"DNA_{dna_id}_All.txt"

    df = clean_data(df)

    json_data = generate_sequences(df)

    G = construct_graph(json_data, k)

    plot_graph(G, graph_filename)

    if is_valid_graph(G):
        dna_sequence = construct_dna_sequence(G)
        output_text = dna_sequence
        all_dna_sequences = construct_all_dna(G)
        save_all_dna(all_dna_sequences, output_all_filename)
        if len(all_dna_sequences) > 1:
            with open(output_all_filename, "r") as f:
                sequences = f.read().strip().split("\n")
            for z, dna in enumerate(sequences, start=1):
                alignment_output_filename = f"DNA_{dna_id}_alignment_{z}.txt"
                align_segments(json_data, dna, alignment_output_filename)
        else:
            alignment_output_filename = f"DNA_{dna_id}_alignment_{1}.txt"
            align_segments(json_data, dna_sequence, alignment_output_filename)
    else:
        output_text = "DNA sequence can not be constructed."

    save_output(output_text, output_filename)


if __name__ == "__main__":
    main()
