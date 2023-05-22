"""
Usage:
python clusters.py -M C00183 -P autism
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def retrieve_metab_bacteria_(metabolite, df_metab_clust, df_bact_clust,
                             df_int_clusters):
    """
    This function saves a bacteria set related to the given metabolite.
    :param metabolite: name of the target metabolite
    :param df_metab_clust: dataframe with 2 cols - name of the metabolite and a corresponding cluster
    :param df_bact_clust: dataframe with 2 cols - name of the bacteria and a corresponding cluster
    :param df_int_clusters: dataframe with cols=metabolite clusters and rows=bacteria clusters,
    df_int_clusters[i][j] represents the correlation coefficient between two clusters

    :return: None
    """
    metab_cl = df_metab_clust[df_metab_clust['metab_name'] == metabolite]['cluster'].values[0]
    cor_clusters = [int(x[0]) for x in df_int_clusters.columns[df_int_clusters.loc[f'{metab_cl}_metab'] > 0]]
    bact_list = df_bact_clust[df_bact_clust['cluster'].isin(cor_clusters)]['bact_name'].values

    with open('bacteria_metab.txt', 'w') as file:
        for bacteria in bact_list:
            file.write("%s\n" % bacteria)

    print('You can find your list of metabolic-based bacterial species in bacteria_metab.txt file')
    return


def cluster_matrix_prep_(metab_clusters, bact_clusters, bact_metab_matrix):
    """
    This function prepares clusters interaction matrix and saves a corresponding heatmap to the current folder.
    :param metab_clusters: dataframe metabolite - its cluster matrix
    :param bact_clusters: dataframe bacteria - its cluster matrix
    :param bact_metab_matrix: dataframe, correlation matrix between bacteria and metabolites

    :return: df_cl_cl: dataframe, bacteria cluster - metabolite cluster interaction
    """

    df_cl_cl = pd.DataFrame(np.zeros((len(metab_clusters['cluster'].unique()), len(bact_clusters['cluster'].unique()))),
                            columns=[str(x) + '_bact' for x in bact_clusters['cluster'].unique()],
                            index=[str(x) + '_metab' for x in metab_clusters['cluster'].unique()])

    for cluster_b in bact_clusters['cluster'].unique():
        bact = list(bact_clusters[bact_clusters['cluster'] == cluster_b]['bact_name'].values)
        for cluster_m in metab_clusters['cluster'].unique():
            metab = list(metab_clusters[metab_clusters['cluster'] == cluster_m]['metab_name'].values)

            res = bact_metab_matrix[bact_metab_matrix['bact_name'].isin(bact)][metab]  # .mean(axis=0)
            df_cl_cl.loc[f'{cluster_m}_metab', [f'{cluster_b}_bact']] = np.nanmean(res.to_numpy())

    heatmap = sns.heatmap(df_cl_cl, cmap="PiYG")
    fig = heatmap.get_figure()
    fig.savefig("bact_metab_clusters_hm.png")
    print('You can find the heatmap for clusters in bact_metab_clusters_hm.png')

    return df_cl_cl


def dataprep(bact_cluster_path, metab_cluster_path, interact_matrix_path):
    # reading files and checking proper naming
    try:
        bact_clusters = pd.read_csv(bact_cluster_path, index_col=0)
        if bact_clusters.shape[1] != 2:
            print('Check columns in bacterial clusters file!')
            return
        # bact_clusters.columns = ['bact_name', 'cluster']
    except FileNotFoundError:
        print('Check the path to file with bacterial clusters!')
        return

    try:
        metab_clusters = pd.read_csv(metab_cluster_path, index_col=0)
        if metab_clusters.shape[1] != 2:
            print('Check columns in metabolic clusters file!')
            return
        # metab_clusters.columns = ['metab_name', 'cluster']
    except FileNotFoundError:
        print('Check the path to file with metabolic clusters!')
        return

    try:
        bact_metab_matrix = pd.read_csv(interact_matrix_path, index_col=0)
        if (bact_metab_matrix.shape != (bact_clusters.shape[0], metab_clusters.shape[0] + 1)) or ():
            print('Check that names in interaction matrix are equal to names in metabolic and bacterial clusters!')
            return
        # bact_metab_matrix.columns = ['bact_name'] + list(bact_metab_matrix.columns[1:])
    except FileNotFoundError:
        print('Check the path to file with correlations!')
        return

    return (metab_clusters, bact_clusters, bact_metab_matrix)


def main(metabolite_arg, phenotype): 
    try:
        kegg = pd.read_csv(f'./data_for_clusters/{phenotype}/kegg_{phenotype}.tsv', sep='\t')
    except FileNotFoundError:
        print(f'Check, that you have a file kegg_{phenotype} with KEGG metadata for metabolites!')
        return
    
    try:
        metabolite = list(kegg[kegg['KEGG'] == metabolite_arg]['Compound'])[0]
    except IndexError:
        print('This compoud was not significantly correlated or not explored. Please, consult list of avaliable metabolites or use pathways!')
        return
    print(f'Searching for {metabolite}..')
    
    dataframes = dataprep(f'./data_for_clusters/{phenotype}/microbe_clusters.csv',
                          f'./data_for_clusters/{phenotype}/metabolite_clusters.csv',
                          f'./data_for_clusters/{phenotype}/interaction_score_matrix.csv')
    if not dataframes:
        print('Something went wrong, check your datasets!')
        return
    else:
        metab_clusters, bact_clusters, bact_metab_matrix = dataframes
        
    df_cl_cl = cluster_matrix_prep_(metab_clusters, bact_clusters, bact_metab_matrix)
    retrieve_metab_bacteria_(metabolite, metab_clusters, bact_clusters, df_cl_cl)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script allows to select bacterial species that correlate with the given metabolite'
    )

    parser.add_argument('-M', '--metabolite', type=str, help='Metabolite argument, KEGG Compound ID')
    parser.add_argument('-P', '--phenotype', type=str,
                        help='Phenotype, folder with precomputed data. Options fo far: cystic_fibrosis, soil, IBD')

    args = parser.parse_args()
    main(args.metabolite, args.phenotype)
