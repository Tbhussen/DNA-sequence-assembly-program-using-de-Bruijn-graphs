"""
author: [Tamim Hussein]
studentnumber: [2502572]
"""


import json


def align_segments(json_data, dna, alignment_output_filename):
    # Create a copy of the data values to manipulate during alignment
    data = json.loads(json_data)
    data_values = list(data.values())
    num_segments = len(data)
    # Stores the aligned segments
    aligned_segments = []
    # Stagging the segments that are currently being aligned
    stacked_segments = []
    alignment = ""
    # Tracking spacing "-"
    spacing = 0
    # Find the first segment that matches the start of the DNA sequence
    # Add it to the aligned segments and stack it for further alignment
    # Remove the matched segment from the data values to avoid reusing it
    for segment in data_values:
        if dna.startswith(segment):
            aligned_segments.append(segment)
            stacked_segments.append(segment)
            alignment += segment + "-" * (len(dna) - len(segment)) + "\n"
            data_values.remove(segment)
            break
    while stacked_segments:
        # Stores the segment that matches the end of the previous segment
        # and its index
        temp = ()
        # Remove the last aligned segment from the stack
        # to find the next segment that can be aligned with it
        previous_segment = stacked_segments.pop()
        for i in range(len(previous_segment) - 1):
            stacked_segments = []
            for segment in data_values:
                if (
                    previous_segment.endswith(segment[:i])
                    and segment not in aligned_segments
                ):
                    stacked_segments.append(segment)
            # If there is only one segment that matches, then
            # it is the next segment to align
            # Store it in temp and
            # continue to the next iteration to search for better matches
            # (larger coverage)
            if len(stacked_segments) == 1:
                temp = (stacked_segments[0], i)
        # If there is a match found, align it with the previous segment and
        # add it to the aligned segments
        # Update the spacing and the alignment string accordingly
        if temp:
            aligned_segments.append(temp[0])
            stacked_segments.append(temp[0])
            spacing += len(previous_segment) - temp[1]
            alignment += (
                "-" * spacing
                + temp[0]
                + "-" * (len(dna) - (len(temp[0]) + spacing))
                + "\n"
            )
            data_values.remove(temp[0])
            continue
        # If no match is found, it means that the current segment
        # can not be aligned with any other segment
        if len(stacked_segments) == 0:
            break
    # If the number of aligned segments is less than the number of segments,
    # it means that there are no segments that can be aligned with the sequence
    # Write the non-aligned segments to the output file with a message
    # indicating that there is no possible alignment for the given DNA sequence
    if len(aligned_segments) < num_segments:
        non_aligned = ""
        for segment in data.values():
            if segment not in aligned_segments:
                non_aligned += segment + "\n"
        with open(alignment_output_filename, "w") as f:
            f.write(f"For this DNA sequence\n{dna}\n")
            f.write("there is no possible alignment "
                    "of the following segments\n")
            f.write(non_aligned)
            return
    # Reaching here means that all segments are aligned, so we
    # write the alignment to the output file
    with open(alignment_output_filename, "w") as f:
        alignment += dna + "\n"
        f.write(alignment)
    return


def construct_dna_sequence(graph):
    # find start node (node with out-degree - in-degree = 1)
    start_node = None
    for node in graph.nodes():
        if graph.out_degree(node) - graph.in_degree(node) == 1:
            start_node = node
            break
    if start_node is None:
        start_node = list(graph.nodes())[0]  # any node can be start node
    # Traverse the graph to construct the DNA sequence
    dna_sequence = start_node  # start with the first node
    # Define adjacency list for all nodes
    # while accounting for multiple edges using list for every node
    adj = {node: [] for node in graph.nodes()}
    for L, R in graph.edges():
        adj[L].append(R)
    # Use stack to track the path and visited edges
    stack = [start_node]
    path = []
    # Use stack-based traversal to find Eulerian path
    while stack:
        current_node = stack[-1]
        # If there are outgoing edges, follow one and remove it from the graph
        if adj[current_node]:
            next_node = adj[current_node].pop()
            stack.append(next_node)
        # If no outgoing edges, add to path and backtrack
        else:
            path.append(stack.pop())
    path = path[::-1]  # reverse to get correct order
    for next_node in path[1:]:
        # append the last character of the next node
        dna_sequence += next_node[-1]
    return dna_sequence


def save_all_dna(sequences: set, filename: str) -> None:
    with open(filename, "w") as f:
        for sequence in sorted(sequences):
            f.write(sequence + "\n")


def construct_all_dna(graph) -> set:
    sequences = set()
    total_edges = graph.number_of_edges()

    start_node = None
    for node in graph.nodes():
        if graph.out_degree(node) - graph.in_degree(node) == 1:
            start_node = node
            break
    if start_node is None:
        start_node = list(graph.nodes())[0]  # any node can be start node

    # Store edges with keys for MultiDiGraph
    edges = list(graph.edges(keys=True))

    stack = [(start_node, set(), [start_node])]

    while stack:
        current_node, used_edges, path = stack.pop()

        if len(used_edges) == total_edges:
            dna_sequence = path[0]

            for node in path[1:]:
                dna_sequence += node[-1]

            sequences.add(dna_sequence)
            continue

        for u, v, key in edges:
            edge_id = (u, v, key)

            if u == current_node and edge_id not in used_edges:
                new_used = used_edges.copy()
                new_used.add(edge_id)

                new_path = path.copy()
                new_path.append(v)

                stack.append((v, new_used, new_path))

    return sequences
