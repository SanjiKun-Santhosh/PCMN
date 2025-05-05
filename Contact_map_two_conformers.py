#!/usr/bin/env python
import matplotlib.pyplot as plt
import mdtraj as md
import contact_map as cm
import networkx as nx
import argparse

def find_contacts(trajectoryA, trajectoryB, filenameA, filenameB,cut_off):
    trajA_contacts = cm.ContactFrequency(trajectoryA,cutoff=cut_off)
    trajB_contacts = cm.ContactFrequency(trajectoryB, cutoff=cut_off)
    plot_contact_map(trajA_contacts, filenameA)
    plot_contact_map(trajB_contacts, filenameB)
    diff = cm.AtomMismatchedContactDifference(trajB_contacts, trajA_contacts)
    fig3, ax3 = diff.residue_contacts.plot(figsize=(8, 6), dpi=300)
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    fig3.savefig(f"diff_contacts_{filenameA}_{filenameB}.pdf")
    fig3.savefig(f"diff_contacts_{filenameA}_{filenameB}.eps")
    diff.residue_contacts.most_common()[:10]
    graphA = trajA_contacts.residue_contacts.to_networkx()
    graphB = trajB_contacts.residue_contacts.to_networkx()
    return graphA, graphB

def plot_contact_map(contacts, filename):
    json_str = contacts.to_json()
    from_json_str = cm.ContactFrequency.from_json(json_str)
    dict = from_json_str.to_json()
    with open(f"sample_{filename}.json", "w") as outfile:
        outfile.write(dict)
    figA, axA = contacts.residue_contacts.plot(figsize=(8, 6), dpi=300)
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    figA.savefig(f"{filename}_contacts.pdf")
    figA.savefig(f"{filename}_contacts.eps")
    contacts.save_to_file(f"{filename}_contacts.p")

def extract_subgraph_for_residue(target_residue, graph):
    target_node = target_residue
    for node in graph:
        node_seq = node.resSeq
        node_name = node.name
        node_label = node_name + str(node_seq)
        if target_node == node_label:
            target_node = node
            break
    if target_node in graph:
        new_graph = nx.Graph()
        new_graph.add_node(target_node)
        for node in graph.nodes():
            if graph.has_edge(target_node, node):
                new_graph.add_edge(target_node, node)
                new_graph.add_node(node)
        return new_graph
    return None

def plot_multiple_contact_graphs(graph_dict, title_prefix="Subgraph for", cols=3, fname_prefix="subgraph"):
    num_graphs = len(graph_dict)
    rows = (num_graphs + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
    axes = axes.flatten()
    for i, (res_label, graph) in enumerate(graph_dict.items()):
        ax = axes[i]
        pos = nx.spring_layout(graph)
        nx.draw(graph, pos, ax=ax, with_labels=True,
                node_size=500, node_color='skyblue', font_size=10, edge_color='gray')
        ax.set_title(f"{title_prefix} {res_label}")
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    plt.savefig(f"{fname_prefix}_combined.png")
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform contact_map analysis")
    parser.add_argument("pdbA_file", help="Path to the PDB file 1")
    parser.add_argument("pdbB_file", help="Path to the PDB_file 2")
    parser.add_argument("cut_off",type=float, default=0.35,help="Cutoff for distance. default is 0.35")
    args = parser.parse_args()

    trajA = md.load(args.pdbA_file, top=args.pdbA_file)
    trajB = md.load(args.pdbB_file, top=args.pdbB_file)
    cut_off=args.cut_off
    filenameA = f"{'.'.join(args.pdbA_file.split('.')[:-1])}"
    filenameB = f"{'.'.join(args.pdbB_file.split('.')[:-1])}"
    
    graphA, graphB = find_contacts(trajA, trajB, filenameA, filenameB,cut_off)

    residues_list = ["LYS171", "HIS201", "TYR202", "LEU203", "GLY204", "LYS205", "GLU239", "ASP258", "GLN395", "KHB360"]

    subgraphs_A = {}
    subgraphs_B = {}
    for res in residues_list:
        subgraphA = extract_subgraph_for_residue(res, graphA)
        subgraphB = extract_subgraph_for_residue(res, graphB)
        if subgraphA:
            subgraphs_A[res] = subgraphA
        if subgraphB:
            subgraphs_B[res] = subgraphB

    plot_multiple_contact_graphs(subgraphs_A, fname_prefix=filenameA)
    plot_multiple_contact_graphs(subgraphs_B, fname_prefix=filenameB)
