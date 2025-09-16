# miRama
Scripts for microRNA analysis from NGS (work in progress)

## miRscreen.pl

Map reads against different small rna datasets and provide evidence that a read is actually that a of a miRNA by excluding from results all reads that don't map to mirna only.
Outputs a tab separated file relating each Unique Molecular Identifier (UMI) to the target miRNA and the relative matches to each library type. 

### Usage

```
-b, --main-bam          bam file containing reads to check
-B, --library-bam       bam files to match the reads against separated by ':' with a type flag to record matches
-v, --verbose           enable verbose mode, Default: 'False'
-h, --help              display help
```

Example:

```
perl mirScreen.pl -v -b bam1 -B bam2:trna -B bam3:mature
```
```
perl mirScreen.pl -v -b bam1 -B bam2:trna,bam3:mature
```

## miRseed.pl

Load gene target sequences from a FASTA file (human_3_UTR.fasta, etc.), mature miRNA sequences from a database FASTA (mature.fa) and Filter miRNAs of interest either from a file containing a list of names, or a regex-like organism prefix (e.g., hsa, bta).

miRNA subsequences can be scanned against the gene targets with three different modes (seed, plus, or minus):

- seed: only extract a subsequence by defining a range of nucleotides (e.g. 2,8 defines a seed starting from the second nucleotide and ending to the 8th, inclusively).

- plus: slide through the miRNA in chunks of a defined size and extract matches.

- minus: reverse-complement the miRNA first, then slide through the miRNA in chunks of a defined size and extract matches.

Outputs two tab separated text tables:

- One relates miRNA IDs to the list of target gene symbols predicted

- One relates Gene Symbols to the list of miRNA IDs targeting them.

### Usage 

```
-t, --targets		target_genes.fasta,
-m, --mirnas		microRNA_database.fasta,
-f, --filter		list of miRNA IDs or regex-like organism prefix (e.g., hsa, bta),
-M, --mode			mode selected: seed,plus or minus, Default: seed,
-s, --slice 		slice size for plus and minus modes, Default: 6,
-r, --range 		range for seed mode, Default: 2,8,
-o, --out-folder 	output directory, Default: ./,
-v, --verbose 		enable verbose mode, Default: False
```
Example: 
```
perl mirSeed.pl -v --targets human_3_UTR.fasta --mirnas mature.fa --filter DE_miRNA.names --mode seed --range 2,7
```

> [!NOTE]
>- "human_3_UTR.fasta" refers to a fasta file containing Homo sapiens 3' UTR sequences downloaded from Ensembl with this [query](http://www.ensembl.org/biomart/martview/01b8569dc51308a39fca679babd4a2c8?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.3utr|hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.external_gene_name|hsapiens_gene_ensembl.default.sequences.description&FILTERS=&VISIBLEPANEL=attributepanel) performed on the "Human genes (GRCh38.p14)" with, in order, the following attributes: "Gene stable ID", "Transcript stable ID", "Gene name", "Gene description".
>- "mature.fa" refers to a fasta file containing mature miRNA sequences downloaded from [miRBase](https://www.mirbase.org/download/).
