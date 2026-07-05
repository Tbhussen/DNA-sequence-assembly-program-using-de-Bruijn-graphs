"""
author: [Tamim Hussein]
studentnumber: [2502572]
"""

import json
from dna_constructor import construct_all_dna
from graph_builder import construct_graph
from project import align_segments, clean_data, is_valid_graph
from project import generate_sequences, construct_dna_sequence
from pytest import mark
import pandas as pd
import networkx as nx


@mark.parametrize(
    "dna_df, expected",
    [
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                    [2, 3, 1, 0, 0, 0],
                    [2, 3, 1, 0, 0, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            pd.DataFrame(
                data=[[2, 1, 1, 0, 0, 0], [2, 2, 0, 1, 0, 0],
                      [2, 3, 1, 0, 0, 0]],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 1, 0, 0],
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                    [3, 1, 1, 0, 0, 1],
                    [3, 2, 0, 1, 0, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            pd.DataFrame(
                data=[[1, 1, 1, 0, 0, 0], [1, 2, 0, 1, 0, 0]],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 0, 0, 0, 1],
                    [1, 2, 0, 0, 0, 1],
                    [2, 1, 1, 0, 0, 0],
                    [2, 1, 0, 1, 0, 0],
                    [2, 2, 0, 0, 1, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            pd.DataFrame(
                data=[[1, 1, 0, 0, 0, 1], [1, 2, 0, 0, 0, 1]],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
        ),
    ],
)
def test_clean_data(dna_df: pd.DataFrame, expected: pd.DataFrame) -> None:
    assert clean_data(dna_df).equals(expected)


@mark.parametrize(
    "dna_df, expected_json_str",
    [
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 2, 0, 1, 0, 0],
                    [1, 3, 0, 0, 1, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            '{"1": "ACG"}',
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 4, 0, 1, 0, 0],
                    [1, 3, 0, 0, 1, 0],
                    [1, 2, 0, 0, 0, 1],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            '{"1": "ATGC"}',
        ),
        (
            pd.DataFrame(
                data=[
                    [2, 1, 1, 0, 0, 0],
                    [1, 2, 0, 1, 0, 0],
                    [2, 3, 0, 0, 1, 0],
                    [1, 4, 0, 0, 1, 0],
                    [1, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                    [1, 3, 0, 0, 1, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            '{"1": "ACGG", "2": "ACG"}',
        ),
        (
            pd.DataFrame(
                data=[
                    [1, 1, 1, 0, 0, 0],
                    [1, 3, 0, 1, 0, 0],
                    [1, 2, 0, 0, 1, 0],
                    [1, 4, 0, 0, 1, 0],
                    [2, 1, 1, 0, 0, 0],
                    [2, 2, 0, 1, 0, 0],
                    [2, 3, 0, 0, 0, 1],
                    [2, 4, 0, 0, 1, 0],
                    [2, 5, 0, 0, 1, 0],
                    [2, 6, 0, 0, 1, 0],
                ],
                columns=["SegmentNr", "Position", "A", "C", "G", "T"],
            ),
            '{"1": "AGCG", "2": "ACTGGG"}',
        ),
    ],
)
def test_generate_sequences(
    dna_df: pd.DataFrame, expected_json_str: str
        ) -> None:
    assert json.loads(generate_sequences(dna_df)
                      ) == json.loads(expected_json_str)


@mark.parametrize(
    "json_data, k, expected_edge_list",
    [
        (
            '{"1": "ATTACTC"}',
            5,
            [("ATTA", "TTAC"), ("TTAC", "TACT"), ("TACT", "ACTC")],
        ),
        (
            '{"1": "ATTACTCGCTAA"}',
            5,
            [
                ("ATTA", "TTAC"), ("TTAC", "TACT"), ("TACT", "ACTC"),
                ("ACTC", "CTCG"), ("CTCG", "TCGC"), ("TCGC", "CGCT"),
                ("CGCT", "GCTA"), ("GCTA", "CTAA")
            ],
        ),
    ],
)
def test_construct_graph(
    json_data: str, k: int, expected_edge_list: list
) -> None:

    graph = construct_graph(json_data, k)
    assert list(graph.edges()) == expected_edge_list


@mark.parametrize(
    "DNA_edge_list,  expected_validity",
    [
        (
            [
                ("ATTA", "TTAC"),
                ("TTAC", "TACT"),
                ("TACT", "ACTC"),
                ("ACTC", "ATTA")
            ],
            True,
        ),
        (
            [("ATTA", "TTAC"), ("TTAC", "TACT"), ("GGGG", "GGGA")],
            False,
        ),
        (
            [
                ("ATTA", "TTAC"),
                ("TTAC", "TACT"),
                ("TACT", "ACTC"),
                ("ACTC", "TTAC"),
                ("TTAC", "ATTA"),
                ("ATTA", "GGGG"),
                ("ATTA", "CCCC"),
            ],
            False,
        ),
        (
            [("ATTA", "TTAC"), ("TACT", "ACTC")],
            False,
        ),
        (
            [
                ("ATTA", "TTAC"),
                ("TTAC", "TACT"),
                ("TACT", "ACTC"),
                ("ACTC", "CTCG"),
                ("CTCG", "TCGC"),
                ("TCGC", "CGCT"),
                ("CGCT", "GCTA"),
                ("GCTA", "CTAA"),
            ],
            True,
        ),
    ],
)
def test_is_valid_graph(DNA_edge_list: list, expected_validity: bool) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert is_valid_graph(debruijn_graph) is expected_validity


@mark.parametrize(
    'DNA_edge_list,  possible_dna_sequence',
    [
        (
                [('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAC')],
                ["AAACAC"]
        ),
        (
            [
                ('CTGA', 'TGAA'), ('TGAA', 'GAAT'), ('GAAT', 'AATG'),
                ('AATG', 'ATGA'), ('ATGA', 'TGAA')
            ],
            ["CTGAATGAA"]
        )
    ])
def test_construct_dna_sequence(
        DNA_edge_list: list, possible_dna_sequence) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert construct_dna_sequence(debruijn_graph) in possible_dna_sequence


@mark.parametrize(
    "DNA_edge_list,  possible_dna_sequences",
    [
        (
            [
                ("CTAA", "TAAT"),
                ("TAAT", "AATT"),
                ("AATT", "ATTA"),
                ("ATTA", "TTAC"),
                ("TTAC", "TACT"),
                ("TACT", "ACTC"),
                ("ACTC", "CTCA"),
                ("CTCA", "TCAC"),
                ("TCAC", "CACT"),
                ("TCAC", "CACT"),
                ("CACT", "ACTG"),
                ("CACT", "ACTG"),
                ("CACT", "ACTA"),
                ("ACTG", "CTGG"),
                ("CTGG", "TGGG"),
                ("TGGG", "GGGT"),
                ("GGGT", "GGTC"),
                ("GGTC", "GTCA"),
                ("GTCA", "TCAC"),
                ("ACTA", "CTAC"),
                ("CTAC", "TACG"),
                ("TACG", "ACGC"),
                ("ACGC", "CGCA"),
                ("CGCA", "GCAC"),
                ("GCAC", "CACT"),
            ],
            {"CTAATTACTCACTACGCACTGGGTCACTG",
             "CTAATTACTCACTGGGTCACTACGCACTG"},
        )
    ],
)
def test_construct_all_dna(
    DNA_edge_list: list, possible_dna_sequences
        ) -> None:
    debruijn_graph = nx.MultiDiGraph()
    for edge in DNA_edge_list:
        debruijn_graph.add_edge(edge[0], edge[1])

    assert construct_all_dna(debruijn_graph) == possible_dna_sequences


@mark.parametrize(
    "json_data, dna_sequence, expected_alignment",
    [
        ('{"1": "ATGAA", "2": "CTGAATGA"}', "CTGAATGAA",
         "CTGAATGA-\n----ATGAA\nCTGAATGAA\n"),
        (
            '{"1": "GGGT", "2": "TTGG", "3": "GTTT"}',
            "GGGTTTGG",
            "GGGT----\n--GTTT--\n----TTGG\nGGGTTTGG\n",
        ),
        (
            '{"1": "ATTA", "2": "TTAC", "3": "TACT", "4": "ACTC"}',
            "GGGGGGGGG",
            "For this DNA sequence\nGGGGGGGGG\nthere is no possible "
            "alignment of the following segments\nATTA\nTTAC\nTACT\nACTC\n",
        ),
        (
            '{"1": "GGGT", "2": "TTGG", "3": "GTTT"}',
            "GGTTTGGG",
            "For this DNA sequence\nGGTTTGGG\nthere is no possible alignment "
            "of the following segments\nGGGT\nTTGG\nGTTT\n",
        ),
    ],
)
def test_align_segments(
            json_data: str, dna_sequence: str, expected_alignment: str
        ) -> None:
    alignment_output_filename = "test_alignment_output.txt"
    align_segments(json_data, dna_sequence, alignment_output_filename)

    with open(alignment_output_filename, "r") as f:
        alignment_result = f.read()

    assert alignment_result == expected_alignment
