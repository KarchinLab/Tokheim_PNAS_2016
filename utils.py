import pandas as pd
from scipy import stats
import numpy as np
import yaml

# logging import
import warnings
import os
import sys

def fetch_significant_genes(input_dir, qval, config):
    """Read the significant driver genes for each method."""
    signif_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        # figure out which columns are being used
        if is_valid_config(config, method_name, 'qvalue'):
            qval_cols = config[method_name]['qvalue']
        else:
            qval_cols = ['qvalue']

        # figure out if a custom score or q-value is used
        if is_valid_config(config, method_name, 'threshold'):
            score_val = float(config[method_name]['threshold']['score'])
            tmp_genes = set()
            for qval_col in qval_cols:
                # get the top scoring genes
                top_direction = config[method_name]['threshold']['top']
                if top_direction == 'low':
                    signif_df = df[df[qval_col]<=score_val]
                else:
                    signif_df = df[df[qval_col]>=score_val]
                tmp_genes |= set(signif_df['gene'].tolist())
            signif_dict[method_name] = list(tmp_genes)
        else:
            # use q-value for threshold
            tmp_genes = set()
            for qval_col in qval_cols:
                # get the significant genes
                signif_df = df[df[qval_col]<=qval]
                tmp_genes |= set(signif_df['gene'].tolist())
            signif_dict[method_name] = list(tmp_genes)
    return signif_dict


def read_method_overlap_genes(path, min_methods):
    # read data
    df = pd.read_table(path)

    # add up all genes with min overlap
    gene_list = []
    for i in range(min_methods, len(df.columns)+1):
        tmp = df[str(i)].dropna().tolist()
        gene_list += tmp

    return gene_list


def process_cgc(path):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    s = df['Mutation Types']
    is_small = s.str.contains('Mis|F|N|S').fillna(False)
    is_somatic = ~df['Tumour Types(Somatic)'].isnull()
    df = df[is_small & is_somatic].copy()

    # get gene names
    cgc_genes = df['Gene Symbol'].tolist()

    return cgc_genes


def fetch_raw_dataframes(input_dir):
    data_dict = {}
    for method_file in os.listdir(input_dir):
        method_name = os.path.splitext(method_file)[0]

        # read in data
        full_path = os.path.join(input_dir, method_file)
        df = pd.read_table(full_path)

        data_dict[method_name] = df

    return data_dict


def fetch_filtered_dataframes(input_dir, output_dir, min_methods, cgc_path=None):
    """Return the result files for each method as a dataframe, but with certain
    genes filtered out. This includes those agreed upon by some minimum
    number of methods or from the CGC.
    """
    # get the raw results
    df_dict = fetch_raw_dataframes(input_dir)

    # get the overlapping genes
    overlap_path = os.path.join(output_dir, 'gene_overlap_counts.txt')
    overlap_genes = read_method_overlap_genes(overlap_path, min_methods)

    # add in cgc genes, if available
    if cgc_path is not None:
        cgc_genes = process_cgc(cgc_path)
        overlap_genes = list(set(overlap_genes + cgc_genes))

    # iterate over each method
    for method in df_dict:
        # filter out agreed upon genes
        df_copy = df_dict[method].copy()
        df_copy = df_copy[~df_copy['gene'].isin(overlap_genes)]

        # update dictionary
        df_dict[method] = df_copy

    return df_dict


def load_config(path):
    """Load YAML configuration file."""
    if path is not None:
        with open(path) as handle:
            config = yaml.load(handle)
        return config
    else:
        return None


def is_valid_config(myconfig, method_name, attribute):
    # check if the above attributes are valid
    has_config = myconfig is not None
    if not has_config:
        return False
    has_method = method_name in myconfig
    if not has_method:
        return False
    has_attribute = attribute in myconfig[method_name]
    if not has_attribute:
        return False

    # no problem, return True
    return True
