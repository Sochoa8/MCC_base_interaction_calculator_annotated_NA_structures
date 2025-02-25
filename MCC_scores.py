#!/usr/bin/env python3
"""
Script to compute Matthews Correlation Coefficient (MCC) scores and other performance metrics for RNA structure predictions.
The user is prompted to select RNApdbee output files (the bpseq file for canonical interactions
and the CSV files for non-canonical and stacking interactions).

Files must follow the naming pattern:
    filtered_<structure_id>_(af|pdb)[_-](2D-bpseq|non-canonical|stacking).<extension>
For example:
    filtered_148d_af-2D-bpseq.txt
    filtered_148d_af-non-canonical.csv
    filtered_148d_af-stacking.csv
    filtered_148d_pdb-2D-bpseq.txt
    filtered_148d_pdb-non-canonical.csv
    filtered_148d_pdb-stacking.csv

The script groups files into pairs by structure and by type (af = prediction; pdb = reference),
computes performance metrics (MCC, precision, recall, F1, and accuracy) for each interaction
category (canonical, non-canonical, stacking) as well as a combined set, and writes the results
to an output CSV file.
"""

import os
import re
import csv
import math
import tkinter as tk
from tkinter import filedialog

# Set DEBUG to True to enable debugging output.
DEBUG = True

# --- Helper to Normalize Header Keys ---
def normalize_key(key):
    """
    Normalize header keys by stripping whitespace, converting to lowercase,
    and replacing various dash-like characters and spaces with a standard hyphen.
    """
    key = key.strip().lower()
    for dash in ["–", "—", "‐"]:
        key = key.replace(dash, "-")
    key = key.replace(" ", "-")
    return key

# ---------------------------
# Parsing Functions (RNApdbee outputs)
# ---------------------------
def parse_canonical(filename):
    """
    Parses a canonical base pairs file in bpseq format.
    Each non-empty, non-comment line is expected to contain at least:
       index  nucleotide  partner_index
    Returns a set of nucleotide pairs (tuples, order normalized).
    """
    interactions = set()
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    i = int(parts[0])
                    j = int(parts[2])
                except ValueError:
                    continue
                if j == 0:
                    continue
                pair = tuple(sorted((i, j)))
                interactions.add(pair)
    except FileNotFoundError:
        print(f"Warning: File not found: {filename}")
    return interactions

def parse_noncanonical(filename):
    """
    Parses a non-canonical interactions CSV file.
    The file is assumed to be semicolon-delimited and to have a header.
    Header keys are normalized using normalize_key so that a column like “Base‐pair”
    becomes "base-pair". Each row's "base-pair" value is expected to be in a format
    like "A.g1 - A.g6". This function extracts the numeric indices and returns a set of pairs.
    """
    interactions = set()
    try:
        with open(filename, 'r', newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            if DEBUG:
                print(f"[DEBUG] Non-canonical raw fieldnames in {os.path.basename(filename)}: {reader.fieldnames}")
            for row in reader:
                row_normalized = {normalize_key(k): (v.strip() if v is not None else "") for k, v in row.items()}
                if DEBUG:
                    print(f"[DEBUG] Normalized row: {row_normalized}")
                bp = row_normalized.get('base-pair')
                if bp:
                    parts = bp.split('-')
                    if len(parts) != 2:
                        continue
                    part1 = parts[0].strip()
                    part2 = parts[1].strip()
                    m1 = re.search(r'\d+', part1)
                    m2 = re.search(r'\d+', part2)
                    if not m1 or not m2:
                        continue
                    try:
                        i = int(m1.group())
                        j = int(m2.group())
                    except ValueError:
                        continue
                    pair = tuple(sorted((i, j)))
                    interactions.add(pair)
        if DEBUG:
            print(f"[DEBUG] Parsed {len(interactions)} non-canonical interactions from {os.path.basename(filename)}")
    except FileNotFoundError:
        print(f"Warning: File not found: {filename}")
    return interactions

def parse_stacking(filename):
    """
    Parses a stacking interactions CSV file.
    Assumes a semicolon-delimited file with a header; keys are normalized.
    Expects a "base-pair" column with values like "A.g1 - A.g6".
    Returns a set of stacking pairs as integer tuples.
    """
    interactions = set()
    try:
        with open(filename, 'r', newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            if DEBUG:
                print(f"[DEBUG] Stacking raw fieldnames in {os.path.basename(filename)}: {reader.fieldnames}")
            for row in reader:
                row_normalized = {normalize_key(k): (v.strip() if v is not None else "") for k, v in row.items()}
                bp = row_normalized.get('base-pair')
                if bp:
                    parts = bp.split('-')
                    if len(parts) != 2:
                        continue
                    part1 = parts[0].strip()
                    part2 = parts[1].strip()
                    m1 = re.search(r'\d+', part1)
                    m2 = re.search(r'\d+', part2)
                    if not m1 or not m2:
                        continue
                    try:
                        i = int(m1.group())
                        j = int(m2.group())
                    except ValueError:
                        continue
                    pair = tuple(sorted((i, j)))
                    interactions.add(pair)
        if DEBUG:
            print(f"[DEBUG] Parsed {len(interactions)} stacking interactions from {os.path.basename(filename)}")
    except FileNotFoundError:
        print(f"Warning: File not found: {filename}")
    return interactions

# ---------------------------
# Additional Metric Functions
# ---------------------------
def compute_precision(TP, FP):
    return TP / (TP + FP) if (TP + FP) > 0 else 0

def compute_recall(TP, FN):
    return TP / (TP + FN) if (TP + FN) > 0 else 0

def compute_f1(precision, recall):
    return 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

def compute_accuracy(TP, TN, FP, FN):
    total = TP + TN + FP + FN
    return (TP + TN) / total if total > 0 else 0

# ---------------------------
# Helper Functions for MCC Calculation
# ---------------------------
def get_nucleotides_count(canonical_file):
    """
    Estimates the number of nucleotides from a canonical file.
    Each non-empty, non-comment line is assumed to correspond to one nucleotide.
    """
    count = 0
    try:
        with open(canonical_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    count += 1
    except FileNotFoundError:
        print(f"Warning: File not found for nucleotide count: {canonical_file}")
    return count

def total_possible_pairs(n):
    """Returns the total number of possible unique pairs for n residues."""
    return (n * (n - 1)) // 2

def compute_counts(ref_set, pred_set, total_possible):
    """
    Computes True Positives (TP), False Positives (FP), False Negatives (FN)
    and True Negatives (TN) given the reference and predicted interaction sets.
    """
    TP = len(ref_set & pred_set)
    FP = len(pred_set - ref_set)
    FN = len(ref_set - pred_set)
    TN = total_possible - (TP + FP + FN)
    return TP, TN, FP, FN

def compute_mcc(TP, TN, FP, FN):
    """Computes the Matthews Correlation Coefficient (MCC)."""
    numerator = (TP * TN) - (FP * FN)
    denominator = math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    if denominator == 0:
        return 0
    return numerator / denominator

# ---------------------------
# File Selection and Grouping
# ---------------------------
def select_and_group_files():
    """
    Opens a file selection dialog so the user can select RNApdbee output files
    (the bpseq file for canonical interactions and the CSV files for non-canonical and stacking).
    Files must follow the naming pattern:
       filtered_<structure_id>_(af|pdb)[_-](2D-bpseq|non-canonical|stacking).<extension>
    Returns a dictionary structured as:
       { structure_id: { 'af': { 'canonical': <path>, 'noncanonical': <path>, 'stacking': <path> },
                         'pdb': { 'canonical': <path>, 'noncanonical': <path>, 'stacking': <path> } } }
    """
    root = tk.Tk()
    root.withdraw()
    filepaths = filedialog.askopenfilenames(
        title="Select RNApdbee output files (bpseq.txt and CSV files)",
        filetypes=[("Text and CSV files", "*.txt *.csv"), ("All Files", "*.*")]
    )
    if not filepaths:
        print("No files selected. Exiting.")
        exit(0)

    files_by_structure = {}
    pattern = re.compile(
        r"filtered_(?P<struct_id>[^_]+)_(?P<type>af|pdb)[_-](?P<category>2D-bpseq|non-canonical|stacking)\.(?P<ext>txt|csv)$",
        re.IGNORECASE
    )
    for fp in filepaths:
        basename = os.path.basename(fp)
        match = pattern.match(basename)
        if not match:
            print(f"File {basename} does not match expected naming pattern. Skipping.")
            continue
        struct_id = match.group("struct_id")
        file_type = match.group("type").lower()
        category = match.group("category").lower()
        if category == "2d-bpseq":
            category = "canonical"
        elif category == "non-canonical":
            category = "noncanonical"
        if struct_id not in files_by_structure:
            files_by_structure[struct_id] = {}
        if file_type not in files_by_structure[struct_id]:
            files_by_structure[struct_id][file_type] = {"canonical": None, "noncanonical": None, "stacking": None}
        files_by_structure[struct_id][file_type][category] = fp

    return files_by_structure

# ---------------------------
# Main Processing Function
# ---------------------------
def process_structure_pairs(structure_files, output_csv):
    """
    For each structure, parses the RNApdbee output files, computes interaction sets,
    calculates MCC, precision, recall, F1, and accuracy by comparing predicted interactions (af) to the reference (pdb),
    and writes the results to the output CSV.
    """
    results = []
    for struct_id, types in structure_files.items():
        print(f"\nProcessing structure: {struct_id}")
        if "af" not in types or "pdb" not in types:
            print(f"  Incomplete pair for structure {struct_id} (missing af or pdb). Skipping.")
            continue

        pred_files = types["af"]
        ref_files = types["pdb"]

        if not pred_files.get("canonical") or not ref_files.get("canonical"):
            print(f"  Missing canonical file for structure {struct_id}. Skipping.")
            continue

        # Parse canonical interactions.
        pred_canonical = parse_canonical(pred_files["canonical"])
        ref_canonical  = parse_canonical(ref_files["canonical"])

        # Parse non-canonical interactions.
        pred_noncanonical = parse_noncanonical(pred_files["noncanonical"]) if pred_files.get("noncanonical") else set()
        ref_noncanonical  = parse_noncanonical(ref_files["noncanonical"]) if ref_files.get("noncanonical") else set()

        # Parse stacking interactions.
        pred_stacking = parse_stacking(pred_files["stacking"]) if pred_files.get("stacking") else set()
        ref_stacking  = parse_stacking(ref_files["stacking"]) if ref_files.get("stacking") else set()

        # --- Debugging Output ---
        if DEBUG:
            print("\n[DEBUG] Non-canonical interactions:")
            print("  Reference non-canonical interactions:")
            print(sorted(ref_noncanonical))
            print("  Predicted non-canonical interactions:")
            print(sorted(pred_noncanonical))
            print("  Intersection (non-canonical):")
            print(sorted(ref_noncanonical & pred_noncanonical))
            print("  Difference (ref - pred, non-canonical):")
            print(sorted(ref_noncanonical - pred_noncanonical))
            print("  Difference (pred - ref, non-canonical):")
            print(sorted(pred_noncanonical - ref_noncanonical))
            
            print("\n[DEBUG] Stacking interactions:")
            print("  Reference stacking interactions:")
            print(sorted(ref_stacking))
            print("  Predicted stacking interactions:")
            print(sorted(pred_stacking))
            print("  Intersection (stacking):")
            print(sorted(ref_stacking & pred_stacking))
            print("  Difference (ref - pred, stacking):")
            print(sorted(ref_stacking - pred_stacking))
            print("  Difference (pred - ref, stacking):")
            print(sorted(pred_stacking - ref_stacking))
            
        # Combine interactions.
        pred_all = pred_canonical | pred_noncanonical | pred_stacking
        ref_all  = ref_canonical  | ref_noncanonical  | ref_stacking

        # Nucleotide count from reference canonical file.
        n_residues = get_nucleotides_count(ref_files["canonical"])
        if n_residues <= 0:
            print(f"  Could not determine nucleotide count for structure {struct_id}. Skipping.")
            continue
        tot_possible = total_possible_pairs(n_residues)

        # --- Compute Metrics for Each Category ---
        # Canonical:
        TP_c, TN_c, FP_c, FN_c = compute_counts(ref_canonical, pred_canonical, tot_possible)
        mcc_canonical = compute_mcc(TP_c, TN_c, FP_c, FN_c)
        precision_canonical = compute_precision(TP_c, FP_c)
        recall_canonical = compute_recall(TP_c, FN_c)
        f1_canonical = compute_f1(precision_canonical, recall_canonical)
        accuracy_canonical = compute_accuracy(TP_c, TN_c, FP_c, FN_c)

        # Non-canonical:
        TP_nc, TN_nc, FP_nc, FN_nc = compute_counts(ref_noncanonical, pred_noncanonical, tot_possible)
        mcc_noncanonical = compute_mcc(TP_nc, TN_nc, FP_nc, FN_nc)
        precision_noncanonical = compute_precision(TP_nc, FP_nc)
        recall_noncanonical = compute_recall(TP_nc, FN_nc)
        f1_noncanonical = compute_f1(precision_noncanonical, recall_noncanonical)
        accuracy_noncanonical = compute_accuracy(TP_nc, TN_nc, FP_nc, FN_nc)

        # Stacking:
        TP_st, TN_st, FP_st, FN_st = compute_counts(ref_stacking, pred_stacking, tot_possible)
        mcc_stacking = compute_mcc(TP_st, TN_st, FP_st, FN_st)
        precision_stacking = compute_precision(TP_st, FP_st)
        recall_stacking = compute_recall(TP_st, FN_st)
        f1_stacking = compute_f1(precision_stacking, recall_stacking)
        accuracy_stacking = compute_accuracy(TP_st, TN_st, FP_st, FN_st)

        # Combined:
        TP_all, TN_all, FP_all, FN_all = compute_counts(ref_all, pred_all, tot_possible)
        mcc_combined = compute_mcc(TP_all, TN_all, FP_all, FN_all)
        precision_combined = compute_precision(TP_all, FP_all)
        recall_combined = compute_recall(TP_all, FN_all)
        f1_combined = compute_f1(precision_combined, recall_combined)
        accuracy_combined = compute_accuracy(TP_all, TN_all, FP_all, FN_all)

        results.append({
            'structure': struct_id,
            'n_residues': n_residues,
            # Canonical metrics
            'MCC_canonical': mcc_canonical,
            'Precision_canonical': precision_canonical,
            'Recall_canonical': recall_canonical,
            'F1_canonical': f1_canonical,
            'Accuracy_canonical': accuracy_canonical,
            # Non-canonical metrics
            'MCC_noncanonical': mcc_noncanonical,
            'Precision_noncanonical': precision_noncanonical,
            'Recall_noncanonical': recall_noncanonical,
            'F1_noncanonical': f1_noncanonical,
            'Accuracy_noncanonical': accuracy_noncanonical,
            # Stacking metrics
            'MCC_stacking': mcc_stacking,
            'Precision_stacking': precision_stacking,
            'Recall_stacking': recall_stacking,
            'F1_stacking': f1_stacking,
            'Accuracy_stacking': accuracy_stacking,
            # Combined metrics
            'MCC_combined': mcc_combined,
            'Precision_combined': precision_combined,
            'Recall_combined': recall_combined,
            'F1_combined': f1_combined,
            'Accuracy_combined': accuracy_combined
        })

        if DEBUG:
            print(f"\n[DEBUG] Interaction counts for structure {struct_id}:")
            print(f"  Canonical: ref={len(ref_canonical)} vs pred={len(pred_canonical)}")
            print(f"  Non-canonical: ref={len(ref_noncanonical)} vs pred={len(pred_noncanonical)}")
            print(f"  Stacking: ref={len(ref_stacking)} vs pred={len(pred_stacking)}")
            print(f"  Total interactions (ref): {len(ref_all)}; (pred): {len(pred_all)}")
            print(f"  Total possible pairs (n={n_residues}): {tot_possible}")
    
    if results:
        with open(output_csv, 'w', newline='') as csvfile:
            fieldnames = [
                'structure', 'n_residues',
                'MCC_canonical', 'Precision_canonical', 'Recall_canonical', 'F1_canonical', 'Accuracy_canonical',
                'MCC_noncanonical', 'Precision_noncanonical', 'Recall_noncanonical', 'F1_noncanonical', 'Accuracy_noncanonical',
                'MCC_stacking', 'Precision_stacking', 'Recall_stacking', 'F1_stacking', 'Accuracy_stacking',
                'MCC_combined', 'Precision_combined', 'Recall_combined', 'F1_combined', 'Accuracy_combined'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r)
        print(f"\nResults written to {output_csv}")
    else:
        print("No valid structure pairs processed; no output written.")

# ---------------------------
# Main Function
# ---------------------------
def main():
    print("Select your RNApdbee output files (bpseq.txt and CSV files).")
    structure_files = select_and_group_files()
    if not structure_files:
        print("No valid files found. Exiting.")
        return

    print("\nGrouped structure files:")
    for struct_id, types in structure_files.items():
        print(f"Structure {struct_id}:")
        for t in types:
            print(f"  {t}: {types[t]}")
    
    root = tk.Tk()
    root.withdraw()
    output_csv = filedialog.asksaveasfilename(
        title="Save MCC results as CSV",
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv")]
    )
    if not output_csv:
        print("No output file chosen. Exiting.")
        return

    process_structure_pairs(structure_files, output_csv)

if __name__ == '__main__':
    main()