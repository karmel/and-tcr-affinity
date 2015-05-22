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