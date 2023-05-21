# SAMOVAR: <a href=""><img src="img/Logo_SAMOVAR_pprpl.png" align="right" width="150" ></a>  <h3> Simulator of Artificial Metanomes: Organisms and Viruses from Abundances to Reads </h3> 

*Metagenomics* – key approach for biological community analysis. Many new tools appear regularly and their validation becomes the crucial challenge. 


Here we come up with an artificial data generation tool `SAMOVAR` that aims aims to improve algorithms development and accelerate scientific discoveries.
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
https://github.com/nvaulin/Metagenomes_Generations_Tool.git
cd Metagenomes_Generations_Tool
```

Create a `conda/mamba` environment with necessary packages and activate it:

```bash
conda env create -f environment.yml
conda activate metageneration
```

### Usage

To run the script, just call it from the directory where the tools is located:

```bash
python Metagenome_generation.py -p [PHENOTYPE] ...
```

To perform the test run use the `2_species` phenotype:
```bash
python Metagenome_generation.py -p 2_species
```

With the real baseline phenotypes its better to select `n` core species with the `ncore` option.
To test the pathways correction use the `example_pathways.txt` file:
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
conda remove --name metageneration --all
rm -rf Metagenomes_Generations_Tool
```

### Citation

If you use these tool, please cite as:
- Vaulin N., Chechenina A., Ivanov. A, Uliantsev V. Metagenomes metageneration metatool. Природа, 2023, 443, 27–28
