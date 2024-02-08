# Calculation of the Percentage Of Conserved Proteins

![](https://img.shields.io/badge/nextflow->=20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-ruby-red)
![](https://img.shields.io/badge/can_use-conda/mamba-yellow.svg)
![](https://img.shields.io/badge/can_use-docker-blue.svg)
![](https://img.shields.io/badge/can_use-singularity-orange.svg)
![](https://img.shields.io/badge/licence-GLP3-lightgrey.svg)

[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer) 

__Update 2023/12: One-vs-All comparisons are now possible in genome and protein input mode. Check `--help` message.__

__Update 2023/10: Now using [Diamond](https://www.nature.com/articles/s41592-021-01101-x) instead of Blastp for protein alignments. Thx [@michoug](https://github.com/michoug) for the Pull Request.__

__Update 2023/05: Re-implementation as a [Nextflow pipeline](nextflow.io). Please feel free to report any [issues](https://github.com/hoelzer/pocp/issues)!__

As input use one amino acid sequence FASTA file per genome such as provided by
[Prokka](https://github.com/tseemann/prokka) or genome FASTA files which will be then annotated via [Prokka](https://github.com/tseemann/prokka). 
The pipeline will then calculate all-vs-all pairwise alignments between all protein sequences and use this
information for POCP calculation following [Qin, Xie _et al_. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24706738). For one-vs-all comparisons see below.

You only need `nextflow` and `conda` or `mamba` or `docker` or `singularity` to run the pipeline. I recommend using `docker`. Then install and run the pipeline:

```bash
# get the pipeline code
nextflow pull hoelzer/pocp 

# check availble release versions and development branches, recommend to use latest release
nextflow info hoelzer/pocp 

# get the help page and define a release version. ATTENTION: use latest version. 
nextflow run hoelzer/pocp -r 2.2.0 --help

# example with genome files as input, all-vs-all comparison, performing a local execution and using Docker
nextflow run hoelzer/pocp -r 2.2.0 --genomes $HOME'/.nextflow/assets/hoelzer/pocp/example/*.fasta' -profile local,docker

# example with protein FASTA files as input (e.g. from Prokka pre-calculated), all-vs-all comparison, performing a SLURM execution and using conda
nextflow run hoelzer/pocp -r 2.2.0 --proteins $HOME'/.nextflow/assets/hoelzer/pocp/example/*.faa' -profile slurm,conda

# example with genome files as input, additional genome file to activate one-vs-all comparison, performing a local execution and using Docker
nextflow run hoelzer/pocp -r 2.2.0 --genomes $HOME'/.nextflow/assets/hoelzer/pocp/example/*.fasta' --genome $HOME/.nextflow/assets/hoelzer/pocp/example/Cav_10DC88.fasta -profile local,docker
```

The final output (`pocp-matrix.tsv`) should look like this (here, the resulting TSV was imported into Numbers on MacOS):

![Example output](example_output.png)

The pipeline identifies orthologous proteins between species using DIAMOND in blastp mode. Please note that the original POCP publication used BLASTP for calculating the alignments. However, DIAMOND is not only faster, which is an advantage when calculating POCP values for larger input data sets, but also achieves the sensitivity of BLASTP ([Buchfink (2021)](https://www.nature.com/articles/s41592-021-01101-x)), especially when using the `--ultra-sensitive` mode, which is activated by default in the pipeline. Another study comparing different alignment programs found that DIAMOND offered the best compromise between speed, sensitivity and quality when a sensitivity option other than the default setting was selected ([Hernández-Salmerón and Moreno-Hagelsieb (2020)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07132-6)). Therefore I decided to use DIAMOND as a more modern solution for the alignment calculation in POCP-nf.

If needed, the following parameters used for filtering the DIAMOND results (blastp mode) can be
adjusted:

```bash
--evalue 1e-5
--seqidentity 0.4
--alnlength 0.5
```

**Attention**: Although the user can customize these parameters, I recommend sticking to the original parameters defined by [Qin, Xie _et al_. 2014](https://www.ncbi.nlm.nih.gov/pubmed/24706738) and otherwise clearly indicating any changed parameter options along with the version of POCP-nf used when sharing POCP results!

Please also note that per default an "all-vs-all" comparison is performed based on the provided FASTA files. However, you can also switch to an "one-vs-all" comparison by additionally providing a single genome FASTA via `--genome` next to the `--genomes` input **or** a single protein multi-FASTA via `--protein` next to the `--proteins` input. In both cases, only "one-vs-all" comparisons will be performed. It is also possible to combine `--genomes` with a target `--protein` FASTA for "one-vs-all" and vice versa. 

If you use the POCP Nextflow pipeline, please cite the original POCP study that introduced the metric and the POCP-nf pipeline:

**[Qin, Qi-Long, _et al_. "A proposed genus boundary for the prokaryotes based on genomic insights." Journal of bacteriology 196.12 (2014): 2210-2215.](https://pubmed.ncbi.nlm.nih.gov/24706738/)**

**[Martin Hölzer. "POCP: An automatic Nextflow pipeline for calculating the percentage of conserved proteins in bacterial taxonomy". SOME JOURNAL. Hopefully in 2024.]()**
