import os
import argparse
import yaml
import subprocess
import warnings
import pandas as pd
import numpy as np
from typing import List
from Bio import Entrez, BiopythonDeprecationWarning
from pysat.examples.hitman import Hitman
from scripts.Genomes_db_update import update_genomes, write_multifasta

GENOMES_DIR = 'genomes'


def parse_args():
    parser = argparse.ArgumentParser(
        usage='Metagenome_generation.py [PHENO]  ...',
        description='''Generate Illumina reads for the  described metagenome in the FILE. '''
    )

    parser.add_argument('-p', '--phenotype', default='2_species', nargs='?',
                        help='the base phenotype for metagenome construction ("Health", "HIV")')
    parser.add_argument('-m', '--metagenome_file', default=None, nargs='?',
                        help='read metagenome composition from the file (tsv with species and abundances)')
    parser.add_argument('--pathways', default=None, nargs='?',
                        help='read matebolic pathways to account from the file (each pathway on the new line')
    parser.add_argument('--metabolites', default=None, nargs='?',
                        help='read metabolites, format: KEGG Compound ID (e.g. C07274)')
    parser.add_argument('--c_model', default='primitive', nargs='?',
                        help='model for core metagenome selection ("primitive", "random", "weighted", "shannon")')
    parser.add_argument('--a_model', default='mean', nargs='?',
                        help='model for species abundances selection ("mean", "normal", "lognormal")')
    parser.add_argument('-c', '--n_core', default=None, nargs='?',
                        help='number of core species to leave in metagenome')
    parser.add_argument('-t', '--threads', default=1, help='number of threads (cores)')
    parser.add_argument('-n', '--n_samples', default=1, nargs='?',
                        help='number of generated metagenome samples')
    parser.add_argument('-o', '--out_dir', default='results', nargs='?',
                        help='path to directory to save generated files')
    parser.add_argument('--email', default='example@email.com', nargs='?',
                        help='Email address for Entrez requests')
    parser.add_argument('--api_key', default=None, nargs='?',
                        help='NCBI API key for Entrez requests (if any)')

    return parser.parse_args()


def generate_core_metagenome(total_metagenome: pd.DataFrame, metagenome_size: int, c_model: str,
                             a_model: str) -> pd.DataFrame:
    """
    Generate a core metagenome by selecting the most abundant species from the given genomes.

    Args:
        total_metagenome: A pandas DataFrame with three columns, the first containing taxid, second - species name and third containing abundance levels.
        metagenome_size: The number of species to select for the core metagenome.
        c_model: The model to use for core metagenome selection. Can be one of 'primitive', 'random', 'weighted', 'shannon'.
        a_model: The model to use for species abundances selection. Can be one of 'mean', 'normal', 'lognormal'.
    Returns:
        A pandas DataFrame with the core metagenome.
    """

    if c_model == 'primitive':
        core = total_metagenome.sort_values(by='mean_abundance', ascending=False).head(metagenome_size)
    elif c_model == 'random':
        core = total_metagenome.sample(n=metagenome_size)
    elif c_model == 'weighted':
        core = total_metagenome.sample(n=metagenome_size, weights='mean_abundance')
    elif c_model == 'random_weighted':
        pass
        # TO-DO random weight from normal distr
        # core = total_metagenome.sample(n=metagenome_size, weights=np.random.normal(core['mean_abundance'], core['sd_abundance']))
    elif c_model == 'shannon':
        total_metagenome['shannon_index'] = - total_metagenome['mean_abundance'] * np.log(
            total_metagenome['mean_abundance'])
        core = total_metagenome.sample(n=metagenome_size, weights='shannon_index')
    else:
        raise ValueError(f"Unknown model for core metagenome selection: {c_model}")

    if a_model == 'mean':
        core['abundance'] = core['mean_abundance']
    elif a_model == 'exponential':
        pass  # TO-DO
    elif a_model == 'normal':
        core['abundance'] = np.random.normal(core['mean_abundance'], core['sd_abundance'])
        while (core['abundance'] <= 0).any():
            core.loc[core['abundance'] <= 0, 'abundance'] = np.random.normal(core['mean_abundance'],
                                                                             core['sd_abundance'])
    elif a_model == 'lognormal':
        core['abundance'] = np.random.lognormal(mean=np.log(core['mean_abundance']), sigma=core['sd_abundance'])
    else:
        raise ValueError(f"Unknown model for species abundances selection: {a_model}")

    return core


def do_hits(metagenome: pd.DataFrame, metabolic_needs: list[list], total_metagenome: pd.DataFrame) -> list:
    """
    This function uses the Hitman package to calculate minimal refill required for metabolic pathways specified
    in a given metabolic database.

    Args:
    metagenome (pd.DataFrame): Baseline metagenome (total or core)
    metabolic_needs (List[List[str]]): A list of lists of metabolites required for each metabolic pathway.
    total_metagenome (pd.DataFrame): Total baseline metagenome

    Returns:
    List[str]: A list representing the species needed to account for the specified metabolites.
    """
    metagenome_species = set(metagenome.species.to_list())
    needs_to_hit = []
    for metabolic_need in metabolic_needs:
        metabolic_need = set(metabolic_need)
        if not set.intersection(metabolic_need, metagenome_species):
            needs_to_hit.append(metabolic_need)
    possible_hits = []
    with Hitman(bootstrap_with=needs_to_hit, htype='sorted') as hitman:
        for hs in hitman.enumerate():
            possible_hits.append(hs)
    hits_scores = []
    for hit in possible_hits:
        hits_scores.append(total_metagenome[total_metagenome['species'].isin(hit)]['abundance'].sum() / len(hit))
    max_index = max(range(len(hits_scores)), key=lambda i: hits_scores[i])
    best_hit = possible_hits[max_index]
    return best_hit


def find_minimal_refill(metagenome: pd.DataFrame, metabolites_specified: List[str],
                        pathways_db: pd.DataFrame, total_metagenome) -> List[str]:
    """
    Given a metagenome composition, a list of metabolites, and a pathways database,
    find the minimal set of additional species that are needed to account for the given metabolites.

    Args:
        metagenome (pd.Dataframe): Baseline metagenome (total or core)
        metabolites_specified (List[str]): A list of metabolite names that need to be accounted for.
        pathways_db (pd.DataFrame): A pandas DataFrame containing the pathways database,
            where rows represent metabolites and columns represent pathways.
            The values are boolean indicating whether a metabolite is present in a pathway.
    Returns:
        List[str]: A list representing the species needed to account for the specified metabolites.
    """
    cols = pathways_db.columns

    missing_pathways = set(metabolites_specified) - set(pathways_db.index)
    if missing_pathways:
        raise KeyError(f"The following elements are missing: {list(missing_pathways)}")
    selected_pathways = pathways_db.loc[metabolites_specified].astype('bool')
    metabolic_needs = selected_pathways.apply(lambda x: list(cols[x.values]), axis=1).to_list()
    return do_hits(metagenome, metabolic_needs, total_metagenome)


def append_species_refill(abudances: pd.DataFrame, to_refill: set) -> pd.DataFrame:
    """
    Append species to an existing dataframe of abundances and adjust abundance levels to maintain normalization.

    Args:
        abudances: A pandas DataFrame with three columns, the first containing taxid, second - species name and third containing abundance levels.
        to_refill: A set of species ids and names to add to the abundance dataframe.
    Returns:
        A pandas DataFrame with the new species added and abundance levels adjusted to maintain normalization.
    Raises:
        ValueError: If the input dataframe does not have the expected two columns, or if the second column
            does not contain numeric data.
    """
    abundance_level = abudances.abundance.mean()
    abundances_refill = pd.DataFrame(to_refill)
    abundances_refill['abundance'] = [abundance_level] * len(to_refill)
    abundances_refill.columns = abundances.columns
    abudances_new = pd.concat([abudances, abundances_refill])
    abudances_new['abundance'] = abudances_new.abundance / abudances_new.abundance.sum()
    return abudances_new


def read_pathways(pathways_input: str) -> List[str]:
    """Reads metabolic pathways from a file or a comma-separated string.

    Args:
        pathways_input (str): Path to a file or a comma-separated string containing metabolic pathways.
    Returns:
        List of strings: A list of metabolic pathway names.
    Raises:
        ValueError: If the input is not a valid path to a file or a comma-separated string.

    """
    if os.path.isfile(pathways_input):
        with open(pathways_input, 'r') as f:
            pathways = [line.strip() for line in f.readlines()]
    elif ',' in pathways_input:
        pathways = pathways_input.split(',')
    else:
        raise ValueError('Invalid input. Please provide a path to a file or a comma-separated string.')
    return pathways


if __name__ == '__main__':
    pheno = parse_args().phenotype
    metagenome_file = parse_args().metagenome_file
    pathways = parse_args().pathways
    metabolites = parse_args().metabolites
    n_core = int(parse_args().n_core)
    n_samples = int(parse_args().n_samples)
    n_threads = int(parse_args().threads)
    RESULTS_DIR = parse_args().out_dir
    core_selection_model = parse_args().c_model
    abundance_selection_model = parse_args().a_model
    email = parse_args().email
    api_key = parse_args().api_key

    Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key

    for sample in range(1, n_samples + 1):
        path = os.path.join(RESULTS_DIR, f'sample_{sample}')
        os.makedirs(path, exist_ok=True)

    baseline_abundances = pd.read_csv(os.path.join('baseline_phenotypes', pheno + '.csv'), header=None, sep=',')
    baseline_abundances.columns = ['tax_id', 'species', 'mean_abundance', 'sd_abundance']
    baseline_abundances['mean_abundance'] = baseline_abundances['mean_abundance'] / baseline_abundances[
        'mean_abundance'].sum()

    pathways_db = pd.read_csv(os.path.join('Databases', 'MetaCyc_pathways_by_species.csv'), sep=';',
                              index_col='Pathways')

    metagenome_size = min(n_core, len(baseline_abundances)) if n_core else len(baseline_abundances)
    abundances = generate_core_metagenome(baseline_abundances, metagenome_size,
                                          core_selection_model, abundance_selection_model)

    for sample in range(1, n_samples + 1):
        to_refill = []
        dir = os.path.join(RESULTS_DIR, f'sample_{sample}')

        if metabolites is not None:
            for metabolite in metabolites:
                metabol_cmd = f'python clusters.py -M {metabolite} -P {pheno}'
                result = subprocess.run(metabol_cmd.split())

        if pathways is not None:
            print('Reading required pathways...')
            pathways_specified = read_pathways(pathways)
            to_refill = find_minimal_refill(abundances, pathways_specified, pathways_db, baseline_abundances)

        if os.path.isfile('bacteria_metab.txt'):
            with open('bacteria_metab.txt', 'r') as file:
                for tx_sp in file:
                    to_refill.append(tuple(tx_sp.split(",")))
        if to_refill:
            abundances = append_species_refill(abundances, to_refill)

        prepared_abundances = update_genomes(GENOMES_DIR, abundances, os.path.join(RESULTS_DIR, f'sample_{sample}'),
                                             n_threads)
        wr_code = write_multifasta(prepared_abundances, GENOMES_DIR)

        iss_params = {
            '-g': os.path.join(GENOMES_DIR, 'multifasta.fna'),
            '--abundance_file': os.path.join(RESULTS_DIR, f'sample_{sample}', 'abundances_for_iss.txt'),
            '-m': 'miseq',
            '-o': os.path.join(RESULTS_DIR, f'sample_{sample}', 'miseq_reads'),
            '--cpus': n_threads
        }
        with open('iss_params.yml', 'r') as f:
            yaml_iss_params = yaml.safe_load(f)
            iss_params = iss_params | yaml_iss_params

        iss_cmd = ['iss', 'generate'] + [str(item) for pair in iss_params.items() for item in pair]
        result = subprocess.run(iss_cmd)

        if os.path.exists(os.path.join(GENOMES_DIR, 'multifasta.fna')):
            os.remove(os.path.join(GENOMES_DIR, 'multifasta.fna'))
        if result.returncode == 0:
            print(f'\nThe sample_{sample} was successfully generated.')
        else:
            print(f'\nThe sample_{sample} generation completed with errors.')

    print('\n Metagenomes generation finished!')