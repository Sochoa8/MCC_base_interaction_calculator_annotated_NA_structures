#!/usr/bin/env python3
import csv
import os
import re
import math
import tkinter as tk
from tkinter import filedialog, messagebox
from collections import Counter, defaultdict

def normalize_interaction(interaction):
    """
    Normalizes an interaction string so that equivalent interactions 
    (e.g. "H/W cis" and "W/H cis") are treated as identical.
    Assumes the format is "X/Y orientation" (e.g., "W/H cis" or "S/W trans").
    The canonical order is defined as: W < H < S.
    """
    parts = interaction.split()
    if len(parts) != 2:
        return interaction  # If the format is unexpected, return it unchanged.
    edge_pair, orientation = parts
    edges = edge_pair.split('/')
    if len(edges) != 2:
        return interaction
    # Define ordering: W < H < S
    ordering = {'W': 1, 'H': 2, 'S': 3}
    edge1 = edges[0].upper()
    edge2 = edges[1].upper()
    # Swap edges if needed
    if ordering.get(edge1, 99) <= ordering.get(edge2, 99):
        canonical_pair = f"{edge1}/{edge2}"
    else:
        canonical_pair = f"{edge2}/{edge1}"
    # Normalize orientation to lowercase.
    orientation = orientation.lower()
    return f"{canonical_pair} {orientation}"

def extract_leontis_westhof(csv_file):
    """
    Reads a non-canonical CSV file (semicolon-delimited) and extracts
    the normalized 'Leontis-Westhof' column values.
    Returns a list of interaction strings.
    """
    lw_list = []
    try:
        with open(csv_file, newline='') as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                lw_value = row.get("Leontis-Westhof")
                if lw_value:
                    norm = normalize_interaction(lw_value.strip())
                    lw_list.append(norm)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
    return lw_list

def get_interaction_counts(file_path):
    """Return a Counter of normalized Leontis-Westhof interactions from a CSV file."""
    lw_list = extract_leontis_westhof(file_path)
    return Counter(lw_list)

def calculate_mcc(y_true, y_pred):
    """
    Calculate the Matthews Correlation Coefficient using its definition.
    
    y_true: list of true binary labels (1 if the interaction is present, 0 otherwise)
    y_pred: list of predicted binary labels
    """
    tp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
    tn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 0)
    fp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
    fn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)
    
    numerator = tp * tn - fp * fn
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return 0.0
    return numerator / denominator

def main():
    # Initialize Tkinter.
    root = tk.Tk()
    root.withdraw()

    messagebox.showinfo("Instructions",
                        "Select all non-canonical CSV files at once.\n"
                        "Each file should be named as:\n"
                        "filtered_<PDB_ID>_<af_or_pdb>-non-canonical.csv\n"
                        "Example: filtered_1ABC_af-non-canonical.csv")
    
    file_paths = filedialog.askopenfilenames(
        title="Select all non-canonical CSV files",
        filetypes=[("CSV Files", "*.csv"), ("CSV Files", "*.CSV"), ("All Files", "*.*")]
    )
    if not file_paths:
        messagebox.showwarning("No Files", "No files selected. Exiting.")
        return

    # Group files by PDB ID using the expected filename pattern.
    pattern = re.compile(
        r"filtered_(?P<pdb_id>[^_]+)_(?P<struct_type>af|pdb)-non-canonical\.csv",
        re.IGNORECASE
    )
    grouped = defaultdict(dict)
    for file in file_paths:
        base = os.path.basename(file)
        match = pattern.search(base)
        if not match:
            print(f"File '{base}' does not match expected pattern. Skipping.")
            continue
        pdb_id = match.group("pdb_id").upper()
        struct_type = match.group("struct_type").lower()
        grouped[pdb_id][struct_type] = file

    # Keep only pairs that have both an 'af' and a 'pdb' file.
    pairs = {pdb_id: files for pdb_id, files in grouped.items() if 'af' in files and 'pdb' in files}
    if not pairs:
        messagebox.showwarning("No Pairs", "No paired files found. Exiting.")
        return

    # Process each paired structure to get interaction counts.
    results = {}  # key: pdb_id, value: {'af': Counter, 'pdb': Counter}
    for pdb_id, files in pairs.items():
        af_counts = get_interaction_counts(files['af'])
        pdb_counts = get_interaction_counts(files['pdb'])
        results[pdb_id] = {'af': af_counts, 'pdb': pdb_counts}

    # Prepare a CSV output that lists interactions for each structure.
    structure_csv = "per_structure_interactions.csv"
    with open(structure_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB_ID", "Source", "Interaction", "Count"])
        for pdb_id, counts in sorted(results.items()):
            # Write PDB interactions.
            for interaction, count in counts['pdb'].items():
                writer.writerow([pdb_id, "PDB", interaction, count])
            # Write Alphafold interactions.
            for interaction, count in counts['af'].items():
                writer.writerow([pdb_id, "AlphaFold", interaction, count])
    print(f"\nDetailed per-structure interaction counts written to '{structure_csv}'.")

    # Determine the union of all normalized interaction types observed.
    all_interactions = set()
    for counts in results.values():
        all_interactions.update(counts['af'].keys())
        all_interactions.update(counts['pdb'].keys())
    all_interactions = sorted(all_interactions)

    # For each interaction type, build binary presence vectors across pairs
    # and record observation counts.
    mcc_scores = {}
    observation_details = {}  # key: interaction, value: (obs_pdb, obs_af, total_pairs)
    for interaction in all_interactions:
        y_true = []  # 1 if the interaction is present in the PDB file.
        y_pred = []  # 1 if the interaction is present in the Alphafold file.
        for pdb_id, counts in results.items():
            pdb_presence = 1 if counts['pdb'][interaction] > 0 else 0
            af_presence  = 1 if counts['af'][interaction] > 0 else 0
            y_true.append(pdb_presence)
            y_pred.append(af_presence)
        # Count observations: number of pairs with the interaction in each.
        obs_pdb = sum(y_true)
        obs_af  = sum(y_pred)
        total = len(y_true)
        observation_details[interaction] = (obs_pdb, obs_af, total)
        # For a single pair, decide manually.
        if len(y_true) < 2:
            mcc = 1.0 if y_true[0] == y_pred[0] else 0.0
        else:
            mcc = calculate_mcc(y_true, y_pred)
        mcc_scores[interaction] = mcc

    # Write summary data to a CSV file.
    summary_csv = "interaction_summary.csv"
    with open(summary_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Interaction", "MCC", "Obs_PDB", "Obs_AF", "Total_Pairs"])
        for interaction, score in sorted(mcc_scores.items(), key=lambda x: x[1]):
            obs_pdb, obs_af, total = observation_details[interaction]
            writer.writerow([interaction, f"{score:.3f}", obs_pdb, obs_af, total])
    print(f"Interaction summary written to '{summary_csv}'.")
    
    # Also print a brief summary to the terminal.
    print("\nManual Matthews Correlation Coefficient (MCC) and Observation Details per Leontis-Westhof interaction:")
    for interaction, score in sorted(mcc_scores.items(), key=lambda x: x[1]):
        obs_pdb, obs_af, total = observation_details[interaction]
        print(f"{interaction}: MCC = {score:.3f} | Observations (PDB): {obs_pdb}/{total}, (AF): {obs_af}/{total}")

if __name__ == "__main__":
    main()