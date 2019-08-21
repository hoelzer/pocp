# Calculation of the Percentage Of Conserved Proteins

As input use one amino acid sequence FASTA file per genome such as provided by
[Prokka](https://github.com/tseemann/prokka). Simply annotate all your genomes
with Prokka first and subsequently organize all your FASTA (*.faa) in one folder
(see example data folder in this repository). The script will then calculate all
pairwise alignments between all FASTA files in the provided folder and use this
information for POCP calculation following [Qin, Xie et al.
2014](https://www.ncbi.nlm.nih.gov/pubmed/24706738).

You need ``awk``, ``grep``, and ``blastp`` as installed dependencies (``awk``
and ``grep`` should be already installed on any linux system). 

If you have blastp installed and in your PATH, simply clone this repository and
execute the following test comand:

````
git clone https://github.com/hoelzer/pocp.git
cd pocp
chmod +x pocp.rb

./pocp.rb example/ example/results/ 2
````

_Important_: The first given parameter must be the input dir holding the
``*.faa`` FASTA files, the second parameter must be the output path and the
third parameter must be the number of threads used for blastp searches.

The final output (``results.csv``) should look like this:

![Example output](example_output.png)

If needed, the following parameters used for filtering the blast results can be
adjusted directly in the ruby script:

````
EVALUE = 1e-5
SEQ_IDENTITY = 0.4
ALN_LENGTH = 0.5
````
