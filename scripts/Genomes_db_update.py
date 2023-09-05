import os
import urllib
import gzip
import urllib.request
import pandas as pd
import concurrent.futures
from Bio import Entrez, SeqIO, Seq
from Bio.SeqRecord import SeqRecord


def generate_ncbi_search_terms(tax_id):
    terms_for_search = []
    specie_for_search = f'"txid{tax_id}"' + '[Organism]'
    filter_for_latest = ' AND latest[filter]'
    filter_for_refseq = '" AND latest refseq"[filter]'
    filter_for_genbank = '"AND latest genbank"[filter]'
    filters_for_relevant = [' AND (NOT anomalous[filter])', '']
    filters_for_taxonomy = ['" AND taxonomy check ok"[filter]', '']
    filters_for_complete = ['" AND complete genome"[filter]', '']
    filters_for_db = [filter_for_refseq, filter_for_genbank, filter_for_latest]

    for filter_for_relevant in filters_for_relevant:
        for filter_for_complete in filters_for_complete:
            for filter_for_taxonomy in filters_for_taxonomy:
                for filter_for_db in filters_for_db:
                    terms_temp = [specie_for_search, filter_for_db, filter_for_relevant]
                    terms_for_search.append(''.join(terms_temp) + filter_for_taxonomy + filter_for_complete)
    terms_for_search.append(f'"txid{tax_id}"')
    return terms_for_search


def get_assembly_download_link(terms_for_search, assembly_summary_cols, assembly_status_translation):
    assemblies_summary = pd.DataFrame(columns=assembly_summary_cols)
    for term in terms_for_search:
        assemblies_info = Entrez.read(Entrez.esearch(db="assembly", term=term, retmax=100))
        assembly_ids = assemblies_info['IdList']
        if len(assembly_ids) > 0:
            break
    for assembly_id in assembly_ids:
        assembly_summary = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id, report="full"))
        assembly_summary = pd.DataFrame(assembly_summary['DocumentSummarySet']['DocumentSummary'],
                                        columns=assembly_summary_cols)
        assemblies_summary = pd.concat([assemblies_summary, assembly_summary])
    assemblies_summary.replace({"AssemblyStatus": assembly_status_translation}, inplace=True)
    assemblies_sorted = assemblies_summary.sort_values(['AssemblyStatus', 'ScaffoldN50',
                                                        'Coverage', 'LastUpdateDate', 'ContigN50'],
                                                       ascending=[True] + [False] * 4)
    links = assemblies_sorted[['FtpPath_RefSeq', 'FtpPath_GenBank']].head(1)
    if links.FtpPath_RefSeq.values:
        url = links.FtpPath_RefSeq.values[0]
    elif links.FtpPath_GenBank.values:
        url = links.FtpPath_GenBank.values[0]
    else:
        raise ValueError(f'''No genome assembly found for {specie} (taxid: {tax_id}). 
                Please, download it manually and store it the {genomes_dir} directory as {fna_filename}''')
    label = os.path.basename(url)
    url = url.removeprefix('ftp://')
    url = 'https://' + url + '/' + label + '_genomic.fna.gz'
    return url


def download_genome(url, fna_filename, genomes_dir):
    urllib.request.urlretrieve(url, os.path.join(genomes_dir, fna_filename))


def update_genomes(abundances, genomes_dir, results_dir, n_threads=1):
    if not os.path.isdir(genomes_dir):
        os.makedirs(genomes_dir)
    prepared_metagenome = pd.DataFrame(columns=['specie', 'taxid', 'abundance'])
    assembly_summary_cols = ['RsUid', 'GbUid', 'AssemblyAccession', 'LastMajorReleaseAccession',
                             'LatestAccession', 'ChainId', 'AssemblyName', 'UCSCName', 'EnsemblName',
                             'Taxid', 'Organism', 'SpeciesTaxid', 'SpeciesName', 'AssemblyType',
                             'AssemblyStatus', 'AssemblyStatusSort', 'WGS', 'GB_BioProjects',
                             'GB_Projects', 'RS_BioProjects', 'RS_Projects', 'BioSampleAccn',
                             'BioSampleId', 'Biosource', 'Coverage', 'PartialGenomeRepresentation',
                             'Primary', 'AssemblyDescription', 'ReleaseLevel', 'ReleaseType',
                             'AsmReleaseDate_GenBank', 'AsmReleaseDate_RefSeq', 'SeqReleaseDate',
                             'AsmUpdateDate', 'SubmissionDate', 'LastUpdateDate',
                             'SubmitterOrganization', 'RefSeq_category', 'AnomalousList',
                             'ExclFromRefSeq', 'PropertyList', 'FromType', 'Synonym', 'ContigN50',
                             'ScaffoldN50', 'AnnotRptUrl', 'FtpPath_GenBank', 'FtpPath_RefSeq',
                             'FtpPath_Assembly_rpt', 'FtpPath_Stats_rpt', 'FtpPath_Regions_rpt',
                             'Busco', 'SortOrder', 'Meta']

    assembly_statuses = {'Complete Genome': 1,
                         'Scaffold': 2,
                         'Contig': 3}
    print(f'Checking {len(abundances)} genomes, it may take some time depending on your internet connection')
    genomes_to_download = []
    for entry in range(len(abundances)):
        tax_id, specie, abundance = abundances[['tax_id', 'species', 'abundance']].iloc[entry]
        terms_for_search = generate_ncbi_search_terms(tax_id)
        prepared_metagenome.loc[len(prepared_metagenome)] = [specie, str(tax_id), abundance]
        fna_filename = str(tax_id) + '.fna.gz'
        if fna_filename not in os.listdir(genomes_dir):
            assembly_link = get_assembly_download_link(terms_for_search, assembly_summary_cols, assembly_statuses)
            genomes_to_download.append((assembly_link, fna_filename))

    if genomes_to_download:
        print(f'Downloading {len(genomes_to_download)} genomes')
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [executor.submit(download_genome, url, fna_filename, genomes_dir) for url, fna_filename in
                       genomes_to_download]
            concurrent.futures.wait(futures)

    print(f'Genomes prepared, writing {results_dir}/metagenome_composition.txt file')
    prepared_metagenome.abundance = prepared_metagenome.abundance / prepared_metagenome.abundance.sum()
    prepared_metagenome.to_csv(os.path.join(results_dir, 'metagenome_composition.txt'), sep='\t', index=False,
                               header=False)
    prepared_metagenome[['taxid', 'abundance']].to_csv(os.path.join(results_dir, 'abundances_for_iss.txt'), sep='\t',
                                                       index=False, header=False)
    return prepared_metagenome


def download_mags(mags_links, mags_folder):
    """
    This function takes a list or np.array (column from the dataframe) with links for each MAG.
    It downloads MAGs according to the link column and returns nothing.
    """
    print(f"Downloading MAGs to {mags_folder}...")
    for link in mags_links:
        if link[link.rfind('/'):] not in os.listdir(mags_folder):
            download_genome(link, link[link.rfind('/')+1:], mags_folder)
    print("MAGs are downloaded!")
    return


def write_multifasta(prepared_abudances, genomes_dir):
    print('Combining the multifasta.fna')
    sequences = []
    for tax_id in prepared_abudances.taxid:
        seq_str = ''
        with gzip.open(os.path.join(genomes_dir, f'{tax_id}.fna.gz'), "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str += record.seq.strip('\n')
        seq_record = SeqRecord(seq=seq_str, id=tax_id, name=tax_id, description='')
        sequences.append(seq_record)
    wr_code = SeqIO.write(sequences, os.path.join(genomes_dir, 'multifasta.fna'), "fasta-2line")
    return wr_code


def write_multifasta(prepared_abudances_taxids, genomes_dir, dir, mag_abundances_ids=None, mags_dir=None):
    print('Gathering the multifasta.fna')
    sequences = []

    for tax_id in prepared_abudances_taxids:
        seq_str = ''
        with gzip.open(os.path.join(genomes_dir, f'{tax_id}.fna.gz'), "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str += record.seq.strip('\n')
        seq_record = SeqRecord(seq=seq_str, id=tax_id, name=tax_id, description='')
        sequences.append(seq_record)

    if mag_abundances_ids is not None:
        for mag_id in mag_abundances_ids:
            seq_str=''
            with gzip.open(os.path.join(mags_dir, f'{mag_id}.gff.gz'), "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq_str += record.seq.strip('\n')
            seq_record = SeqRecord(seq=seq_str, id=mag_id, name=mag_id, description='')
            sequences.append(seq_record)

    wr_code = SeqIO.write(sequences, os.path.join(dir, 'multifasta.fna'), "fasta-2line")
    return wr_code
