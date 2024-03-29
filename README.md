# SAMOVAR: <a href=""><img src="img/Logo_SAMOVAR_pprpl.png" align="right" width="150" ></a>  <h3> Simulator of Artificial Metagenomes: Organisms and Viruses from Abundances to Reads </h3> 

*Metagenomics* – key approach for biological community analysis. Many new tools appear regularly and their validation becomes the crucial challenge. 


Here we come up with an artificial data generation tool `SAMOVAR` that aims to improve algorithms development and accelerate scientific discoveries.
This  pipeline takes the input phenotype or environment as a community property and technical parameters for NGS library to generate `fastq` files with Illumina reads. SAMOVAR performes additional bacteria selection  based on a given *metabolite* or *metabolic pathways*. Overall, our pipeline consists of four steps: 
- Core-bacterial set selection by given phenotype;
- Selection of additional bacteria species based on metabolites and metabolic pathways;
- Prediction of the relative abundance for a selected set of bacteria;
- Illumina reads generation. 


<div style='justify-content: center'>
<img src="img/pipeline.png" align='center', width="80%">
</div>

### Installation

To get the tool clone the git repository:

```bash
git clone https://github.com/nvaulin/samovar.git && cd samovar
```

Create a `conda/mamba` environment with necessary packages and activate it:

```bash
conda env create -f environment.yml
conda activate samovar
```

Then update ingredients for the samovar (genomes database, 5.5 Gb):

```bash
wget https://www.dropbox.com/sh/goeh43tyhu62es3/AADT2w0FBB8kT1z0dXqt-UHKa\?dl\=0 -O genomes_files && unzip -o genomes_files -d genomes
```

### Usage

To run the script, just call it from the directory where the tool is located:

```bash
python Metagenome_generation.py -p [PHENOTYPE] ...
```

Usage options:

```bash
options:
  -h, --help            show this help message and exit
  -p [PHENOTYPE], --phenotype [PHENOTYPE]
                        the base phenotype for metagenome construction
                        ("Health", "HIV")
  -m [METAGENOME_FILE], --metagenome_file [METAGENOME_FILE]
                        read metagenome composition from the file (tsv with
                        species and abundances)
  --pathways [PATHWAYS]
                        read matebolic pathways to account from the file (each
                        pathway on the new line
  --metabolites [METABOLITES]
                        read metabolites, format: KEGG Compound ID (e.g.
                        C07274)
  --c_model [C_MODEL]   model for core metagenome selection ("primitive",
                        "random", "weighted", "weighted_lognormal",
                        "weighted_exponential", "shannon")
  --a_model [A_MODEL]   model for species abundances selection ("mean",
                        "exponential", "normal", "lognormal")
  -c [N_CORE], --n_core [N_CORE]
                        number of core species to leave in metagenome
  -t THREADS, --threads THREADS
                        number of threads (cores)
  -n [N_SAMPLES], --n_samples [N_SAMPLES]
                        number of generated metagenome samples
  -r [N_READS], --n_reads [N_READS]
                        number of reads to generate (if set, overwrites the
                        number present in iss_params.yml)
  -o [OUT_DIR], --out_dir [OUT_DIR]
                        path to directory to save generated files
  --email [EMAIL]       Email address for Entrez requests
  --api_key [API_KEY]   NCBI API key for Entrez requests (if any)

```

Additional *InSilicoSeq* reads generation parameters (such as number of reads, error model, etc) can be also specified in the `iss_params.yml` file (the full list of the *InSilicoSeq* parameters can be found in it's [documentation](https://insilicoseq.readthedocs.io/en/latest/)).

### Examples

To perform the test run use the `2_species` phenotype:
```bash
python Metagenome_generation.py -p 2_species
```

With the real baseline phenotypes its better to select *n* core species with the `ncore` (`c`) option. We **HIGHLY recommend always use `ncore` species** when working with real metagenomes. To test the pathways correction use the `example_pathways.txt` file:
```bash
python Metagenome_generation.py -p Health -c 10 --pathways example_pathways.txt
```

To get more information about the particular script, run:

```bash
python Metagenome_generagtion.py  --help
```


### Uninstallation

To uninstall the tool remove the conda environment and delete the cloned folder:
```python
conda remove --name samovar --all
rm -rf samovar
```

### Citation

If you use these tool, please cite as:
- Chechenina A., Vaulin N., Ivanov. A, Ulyantsev V. Development of *in-silico* models of metagenomic communities with given properties and a pipeline for their generation. Bioinformatics institute 2022/23, 2023, 22-24

```bibtex
@article{samovar2023,
  title={Development of in-silico models of metagenomic communities with given properties and a pipeline for their generation},
  author={Chechenina, A. and Vaulin, N. and Ivanov, A. and Ulyantsev, V.},
  journal={Bioinformatics institute 2022/23},
  year={2023},
  pages={22--24}
}
```
