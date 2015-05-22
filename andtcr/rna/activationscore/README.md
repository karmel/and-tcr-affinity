## Using the Activation Signature Score

The activation signature score, as described in the related manuscript, can be used to rank arbitrary CD4+ T cell data sets. Note that it works best for mouse genes at this point, but could be adapted to work for human genes if the names were properly matched up.

The simplest path to use:

1. Download the [calculate_activation_signature_score.py](calculate_activation_signature_score.py) script and the list of genes in [activation_score_genes](activation_score_genes.txt) into a directory.

1. Prepare the input file. For a given set of samples, you should generate a **tab-delimited** file:
	
	- The first column should be labeled either `gene symbol` or `refseq`, and should contain either the standard gene symbol (for example, `Tnf`) for the given transcript or the murine RefSeq identifier (for example, `NM_013693`).

	- Each subsequent column should contain normalized expression values (RPKM or microarray probe intensity) for a given sample, with the sample name in the header for the column.

	- There are many examples of preparing such files in [the examples directory](examples); here is an example of a properly formatted file:

    gene symbol		sample_1	sample_2	sample_3
    Il6ra			10.4		20.2		30.0
    Ninj1			1.6			2.3			0.0
    Iars			12.5		43.6		6.1
    ...

1. Run the script with the input file as a parameter on the command line. The output will be printed to standard out:

    python calculate_activation_signature_score.py -f my_gene_values.txt

    Activation scores:
    naive-1   -0.767299
    naive-2   -0.783955
    th1-1      0.514566
    th1-2      0.548627
    th2-1      0.617183
    th2-2      0.607981
    th17-1     0.346837
    th17-2     0.394038
    itreg-1    0.241824
    itreg-2    0.226667
    ntreg-1   -0.946469
    ntreg-2   -1.000000
    dtype: float64

And those are your activation signature scores to rank samples and plot as you please!

[The examples directory](examples) contains many examples of how to easily create the input files from RPKM or microarray matrices, and also how to use the [calculate_activation_signature_score.py](calculate_activation_signature_score.py) script in Python directly, rather than at the command line, so that the scores can be more easily combined, analyzed, and plotted.

Questions? Comments? Email me: [karmel@arcaio.com](mailto:karmel@arcaio.com)