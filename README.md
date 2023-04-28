# Calculation of the Percentage Of Conserved Proteins

![](https://img.shields.io/badge/uses-ruby-red)
![](https://img.shields.io/badge/licence-GLP3-lightgrey.svg)

[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer) 

__Update 2023: Re-implementation as a [Nextflow pipeline](nextflow.io).__

__Please note that I simply provide this code for POCP calculation and it is not polished in any way. Thus, the user experience might be not great, but the code does the job. Please feel free to report any [issues](https://github.com/hoelzer/pocp/issues)!__

As input use one amino acid sequence FASTA file per genome such as provided by
[Prokka](https://github.com/tseemann/prokka). Simply annotate all your genomes
with Prokka first and subsequently organize all your FASTA (`*.faa`) in one folder
(see example data folder in this repository). The script will then calculate all
pairwise alignments between all FASTA files in the provided folder and use this
information for POCP calculation following [Qin, Xie _et al_.
2014](https://www.ncbi.nlm.nih.gov/pubmed/24706738). You can also compare all your `*.faa`
input files only against one other `*.faa` file to reduce runtime (see below).

You need `ruby`, `awk`, `grep`, and `blastp` as installed dependencies (`ruby`, `awk`
and `grep` should be already installed on any Linux system). The easiest way 
to install `blastp` is via [conda](https://docs.conda.io/en/latest/miniconda.html):

```bash
conda create -n pocp -c bioconda blast && conda activate pocp 
```

If you have blastp installed and in your PATH, simply skip the above step and clone this repository and
execute the following test comand:

```bash
git clone https://github.com/hoelzer/pocp.git
cd pocp
chmod +x pocp.rb

./pocp.rb example/ example/results/ 2
```

_Important_: The first given parameter must be the input dir holding the
``*.faa`` FASTA files, the second parameter must be the output path and the
third parameter must be the number of threads used for `blastp` searches.

The final output (`results.csv`) should look like this:

![Example output](example_output.png)

If needed, the following parameters used for filtering the blast results can be
adjusted directly in the ruby script:

```bash
EVALUE = 1e-5
SEQ_IDENTITY = 0.4
ALN_LENGTH = 0.5
```

Per default all pairwise comparisons of all `.faa` FASTA files located in the input folder are performed. 
Please define a single FASTA filename (not path, only the filename), if comparisons should be only performed against the protein sequences in this file:

```bash
./pocp.rb example/ example/results/ 2 Cav_10DC88.faa
```
