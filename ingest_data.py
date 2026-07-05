"""
author: [Tamim Hussein]
studentnumber: [2502572]
"""

import pandas as pd
import json
import re


def read_csv(name: str) -> tuple[str, int, pd.DataFrame]:
    """
    Reads a CSV file and returns sequence ID, k-mer size,
    and DataFrame.
    """
    pattern = re.search(r'.*DNA_(\d+)_(\d+)', name)
    if pattern:
        seq_id = pattern.group(1)
        seq_k = int(pattern.group(2))
    else:
        raise ValueError("Filename does not match the"
                         " expected pattern 'DNA_<id>_<k>.csv'")
    data = pd.read_csv(name, header=None)
    # Rename columns
    headers = ["SegmentNr", "Position", "A", "C", "G", "T"]
    data.columns = headers
    return seq_id, seq_k, data


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Cleans the DataFrame by removing duplicates and ensuring
    that each position has exactly one nucleotide.
    Solves Errors 1, 2, 3, and 4
    """
    # Remove Duplicates
    df = df.drop_duplicates()
    unique_segments = []
    # Check for unique positions within each segment
    # and valid nucleotide counts
    for segment in df['SegmentNr'].unique():
        unique_nucleotides_per_segment = []
        segment_df = df[df['SegmentNr'] == segment]
        m = max(segment_df['Position'])
        for pos in range(1, m + 1):
            # Remove segments with duplicate or missing positions
            if len(segment_df[segment_df["Position"] == pos]) != 1:
                df = df.drop(df[df["SegmentNr"] == segment].index)
                break
            # Check if exactly one nucleotide is present
            rowSum = sum(
                segment_df[segment_df['Position'] == pos][
                    ['A', 'C', 'G', 'T']
                    ].values[0]
                ) == 1
            # Remove segments with invalid nucleotide counts
            if not rowSum:
                df = df.drop(df[(df['SegmentNr'] == segment)].index)
                break
            unique_nucleotides_per_segment.append(
                segment_df[segment_df["Position"] == pos][
                    ['A', 'C', 'G', 'T']
                    ].values[0].tolist()
                    )
        # Remove segments with duplicate nucleotide patterns
        if unique_nucleotides_per_segment in unique_segments:
            df = df.drop(df[df["SegmentNr"] == segment].index)
        else:
            unique_segments.append(unique_nucleotides_per_segment)
    df = df.reset_index(drop=True)
    return df


def generate_sequences(df: pd.DataFrame) -> str:
    """
    Converts the cleaned DataFrame into a JSON format.
    """
    data = {}
    # Extracts the Sequence from the DataFrame and builds the JSON structure
    for segment in df['SegmentNr'].unique():
        # Sort by position in ascending order to ensure correct
        # sequence construction
        segment_df = df[df['SegmentNr'] == segment].sort_values(by='Position')
        result = ""
        for _, row in segment_df.iterrows():
            result += ('A' if row['A'] == 1
                       else 'C' if row['C'] == 1
                       else 'G' if row['G'] == 1
                       else 'T' if row['T'] == 1
                       else 'N')
        data[str(segment)] = result  # Change type for segment ID
    json_data = json.dumps(data, indent=4)
    return json_data
