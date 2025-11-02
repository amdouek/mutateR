mutateR
================
Alon M Douek
2025-11-03

- [Overview](#overview)
  - [Background](#background)
  - [Why does this matter?](#why-does-this-matter)
  - [How can we prevent this?](#how-can-we-prevent-this)
- [How does `mutateR` work?](#how-does-mutater-work)
- [Getting started](#getting-started)
  - [Installation](#installation)
- [Minimal `mutateR` workflow](#minimal-mutater-workflow)
  - [Interpreting the `run_mutateR()`
    output](#interpreting-the-run_mutater-output)
  - [Special cases](#special-cases)
- [Manual function execution](#manual-function-execution)
  - [`get_gene_info()`](#get_gene_info)
  - [`get_exon_structures()`](#get_exon_structures)
  - [`find_cas9_sites()`](#find_cas9_sites)
  - [`find_cas12a_sites()`](#find_cas12a_sites)
  - [`score_grnas()`](#score_grnas)
- [To be implemented](#to-be-implemented)
- [Session information](#session-information)
- [References](#references)

## Overview

`mutateR` streamlines design of frame-aware CRISPR/Cas-mediated
deletions.

It identifies pairs of non-contiguous exons in a given transcript that
are phase-compatible with each other (*i.e.*, pairs of exons that do not
normally neighbour each other in the transcript, but the 5’ exon ends in
the same reading frame position (0, 1, or 2) as the 3’ exon begins
with).

### Background

The underlying motivation is to facilitate the generation of mutant
alleles that are not prone to nonsense-mediated decay (NMD), which can
induce transcriptional adaptive responses ([El-Brolosy et al., *Nature*
(2019)](https://doi.org/10.1038/s41586-019-1064-z "Genetic compensation triggered by mutant mRNA degradation")
and [Ma et al., Nature
(2019)](https://www.nature.com/articles/s41586-019-1057-y "PTC-bearing mRNA elicits a genetic compensation response via Upf3a and COMPASS components")).
This behaviour can attenuate the phenotypic potency of mutant alleles,
making the study of gene loss-of-function more difficult.

NMD is a complex process that cells use to survey and degrade
transcripts that contain a premature termination codon (PTC). This
generally occurs when a mutation causes a frameshift, resulting in a
termination codon being introduced in a non-terminal exon. During
splicing of mature mRNA, exon junction protein complexes (EJCs) are
deposited at exon-exon junctions; during translation, the ribosome
displaces these complexes. If a ribosome encounters a PTC **at least 50
nucleotides** upstream of the last exon-exon junction, it will fail to
displace any downstream EJCs. The remaining EJC(s) will interact with
termination factors to induce decay of the PTC-bearing transcript.

While the degradation of PTC-bearing mRNAs protects the cell from
producing potentially toxic truncated proteins, the mRNA fragments that
result from NMD can induce upregulation of other (‘adapting’) genes.
This effect is a compensatory behaviour termed *‘transcriptional
adaptation’*, and can cause masking of gene loss-of-function phenotypes.

### Why does this matter?

Suppose you are studying a novel gene of interest (‘A’). Gene A is
expressed in the brain, and you hypothesise that loss of Gene A causes a
brain growth phenotype. When you knock down expression of Gene A in your
model organism of choice, you see the expected stunting of brain
development; however, when you look at the brains of animals homozygous
for a frameshift mutation causing a PTC early in the gene, the brains
look normal. Gene A is functionally related to, and expressed in the
same region as, another gene (‘B’). When you measure expression of Gene
A and Gene B in these brains, you see little to no Gene A expression,
but Gene B is expressed at higher levels than in brains from animals
that do not have this mutation. This is because Gene B has compensated
for the loss of Gene A by responding to the mRNA fragments produced by
NMD of the Gene A transcript. While this is *great* for the animal (it
has a largely “normal” brain!), it means that we can’t really use our
Gene A mutant to reliably study the effects of loss of Gene A function.

### How can we prevent this?

In short, **use tools for targeted mutagenesis to their strengths**.
Back when the only way to make a mutant animal was to randomly
mutagenise its genome and then recover alleles based on phenotypes,
researchers didn’t have a great deal of control over the nature of the
resulting alleles. Modern sequence-targeting methods (such as
CRISPR/Cas9) allow for much more control over the nature of the
mutations that we can generate.

If we delete parts of a gene that correspond to functional domains of
the encoded protein, but maintain the original reading frame (such that
there is no PTC**\***), we are left with a stably-expressed, truncated
mutant transcript that does not encode a functional protein but also
does not induce transcriptional adaptation via NMD.

> *\*Note*: There is nuance to NMD induction. For example, PTCs are much
> better tolerated (*i.e.*, do not induce NMD) if they occur in the
> transcript’s terminal exon. For an excellent quantitative study of NMD
> induction relative to transcript PTC position, please see [Lindeboom,
> Supek and Lehner, *Nature Genetics*
> (2016)](https://doi.org/10.1038/ng.3664 "The rules and impact of nonsense-mediated mRNA decay in human cancers").

This type of deletion can be achieved by simultaneously targeting two
sites in the gene of interest; the intervening region between the two
targeted sites will be excised, and the ends of the genetic lesion are
joined together by DNA repair mechanisms.

## How does `mutateR` work?

`mutateR` uses `biomaRt` to retrieve transcript and exon phase
information for your gene of interest, then calculates which pairs of
non-contiguous exons are phase-compatible. It then scans these exons for
Cas effector PAMs (currently Cas9 NGG and Cas12a TTTV) associated with
protospacers that do not cross exon boundaries. After calculating
on-target scores for potential gRNAs using rule sets from `crisprScore`
([Hoberecht et al., *Nature Communications*
(2022)](https://www.nature.com/articles/s41467-022-34320-7 "A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies")),
`mutateR` returns recommended gRNA pairs to target exons where the
flanking exons are in-frame with each other. It also produces a
graphical representation of the selected transcript, its
phase-compatible exons, and the recommended exon pairs for targeting.

The `mutateR` package consists of the following ordered functions:

| Step | Function | Functionality |
|----|----|----|
| 1 | `get_gene_info()` | Query Ensembl for gene and transcript metadata. |
| 2 | `get_exon_structures()` | Retrieve exon coordinates and phase metadata. |
| 3 | `find_cas9_sites()`/`find_cas12a_sites()` | Identify nuclease recognition sites. |
| 4 | `check_exon_phase()` | Determine phase-compatible exon pairs. |
| 5 | `check_frameshift_ptc()` | Calculate and flag frameshift/PTC consequences. |
| 6 | `map_protein_domains()` | Retrieve protein domain annotations from InterPro corresponding to exons. |
| 7 | `score_gRNAs()` | Score guide RNAs via `crisprScore`. |
| 8 | `filter_valid_gRNAs()` and `assemble_gRNA_pairs()` | Filter allowed gRNA pairs and assemble into a readable data.frame. |
| 9 | `plot_gRNA_design()` | Visualise phase compatibility of non-contiguous exon pairs and top-ranked valid gRNA pairs. |
| 10 | `plot_grna_heatmap()` | For large gene (\>12 exons) visualisation - helper function for `plot_grna_design()` but can be executed separately. |
| 11 | `run_mutateR()` | Wrapper; runs entire pipeline. |

`mutateR` also has an additional function, `design_gRNA_pairs()`, that
is not included in the wrapper function. This function can be used in
isolation to produce a data.frame of gRNA pairs for the gene of
interest.

## Getting started

### Installation

Currently, `mutateR` can be installed from GitHub:

``` r
devtools::install_github("amdouek/mutateR")
library(mutateR)
```

Running `mutateR` requires only installation of the package itself, and
that you have also installed the `BSgenome` R package for your species
of interest (e.g., `BSgenome.Hsapiens.UCSC.hg38` for human). These
BSgenome data packages can be quite large, so are not installed
alongside `mutateR`.

> For convenience’s sake, here are the Bioconductor pages to install the
> genomes for commonly-used species:
>
> - [human
>   (hg38)](https://doi.org/doi:10.18129/B9.bioc.BSgenome.Hsapiens.UCSC.hg38),
>
> - [mouse
>   (mm10)](https://doi.org/doi:10.18129/B9.bioc.BSgenome.Mmusculus.UCSC.mm10),
>
> - [rat
>   (rn7)](https://doi.org/doi:10.18129/B9.bioc.BSgenome.Rnorvegicus.UCSC.rn7),
>
> - [zebrafish
>   (danRer11)](https://doi.org/doi:10.18129/B9.bioc.BSgenome.Drerio.UCSC.danRer11),
>
> - [*Drosophila melanogaster*
>   (dm6)](https://doi.org/doi:10.18129/B9.bioc.BSgenome.Dmelanogaster.UCSC.dm6)
>
> - [*Caenorhabditis elegans*
>   (ce11)](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Celegans.UCSC.ce11.html)
>
> You can also see all genomes currently available through `BSgenome` by
> running `BSgenome::available.genomes()`.

In the below example run, we will run `mutateR` on human *TP53* to
select optimal Cas9 gRNA pairs*.* Assuming you have installed both
`mutateR` and your BSgenome of choice (in this case, human):

``` r
library(mutateR)
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Minimal `mutateR` workflow

The wrapper function `run_mutateR()` executes the full pipeline (using
the canonical transcript by default). For human *TP53:*

``` r
tp53_res <- run_mutateR(
  gene_id = "TP53",
  species = "hsapiens",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  nuclease = "Cas9",
  score_method = 'ruleset1',
  quiet = FALSE,
  plot_mode = 'heat'
  )
#> Retrieving gene/transcript information...
#> Using transcript: ENST00000269305 for gene: TP53
#> Locating Cas9 target sites...
#> Scoring gRNAs using model: ruleset1
#> Computing on‑target scores using ruleset1 model...
#> Re‑orienting 214 guides so PAM = NGG at positions 26–27.
#> PAM triplet distribution (positions 25–27):
#> 
#> AGG CGG GGG TGG 
#> 102  36 136 125
#> Scored 399 guides using ruleset1.
#> Assembling valid gRNA pairs for TP53 ...
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Flattening on-target scores from GRanges ...
#> Retrieving InterPro domain annotations from Ensembl Genes mart...
#> Generated 2664 candidate exon‑flanking gRNA pairs.
#> Plotting exon phase compatibility and gRNA pairs...
#> Retrieving Pfam domain annotations from Ensembl Genes mart...
#> mutateR pipeline completed for TP53, finding 2664 gRNA pairs.
```

> Note that this function also reports a PAM distribution. In the case
> of *TP53*, CGG PAMs are significantly underrepresented amongst those
> identified by `mutateR`, as is generally expected in the human exome
> (see [Yang, Zhu and Jin, *BMC Genomics*
> (2024)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-11073-9)
> and [Scott and Zhang, *Nature Medicine*
> (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5749234/)).

The object `tp53_res` is a list with the following levels:

``` r
names(tp53_res)
#> [1] "gene_id"       "gene_symbol"   "transcript_id" "exons"        
#> [5] "scored_grnas"  "pairs"         "plot"
```

With the actual selected gRNA pairs (with associated metadata) stashed
as a data.frame in `tp53_res$pairs`.

``` r
kableExtra::scroll_box(knitr::kable(head(tp53_res$pairs,5)),width = '700px')
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:700px; ">

<table>

<thead>

<tr>

<th style="text-align:right;">

upstream_pair
</th>

<th style="text-align:right;">

downstream_pair
</th>

<th style="text-align:right;">

exon_5p
</th>

<th style="text-align:right;">

exon_3p
</th>

<th style="text-align:left;">

compatible
</th>

<th style="text-align:left;">

frameshift
</th>

<th style="text-align:left;">

ptc_flag
</th>

<th style="text-align:left;">

terminal_exon_case
</th>

<th style="text-align:left;">

protospacer_sequence_5p
</th>

<th style="text-align:left;">

pam_sequence_5p
</th>

<th style="text-align:right;">

ontarget_score_5p
</th>

<th style="text-align:left;">

protospacer_sequence_3p
</th>

<th style="text-align:left;">

pam_sequence_3p
</th>

<th style="text-align:right;">

ontarget_score_3p
</th>

<th style="text-align:left;">

domains
</th>

<th style="text-align:left;">

recommended
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

ACCATTGTTCAATATCGTCC
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.0442887
</td>

<td style="text-align:left;">

AAACACTTTTCGACATAGTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0507973
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

AGACCTGCCCTGTGCAGCTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0471244
</td>

<td style="text-align:left;">

AAACACTTTTCGACATAGTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0507973
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53-like_TF_DNA-bd_sf; p53_DNA-bd;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

CCTTCCCAGAAAACCTACCA
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5585539
</td>

<td style="text-align:left;">

AAACACTTTTCGACATAGTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0507973
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

CAGGGCAGGTCTTGGCCAGT
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0115891
</td>

<td style="text-align:left;">

AAACACTTTTCGACATAGTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0507973
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53-like_TF_DNA-bd_sf; p53_DNA-bd;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

AAGATGACAGGGGCCAGGAG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.0817923
</td>

<td style="text-align:left;">

AAACACTTTTCGACATAGTG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.0507973
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

FALSE
</td>

</tr>

</tbody>

</table>

</div>

By default, all pairs are present in this data.frame; if you are
interested in only those that `mutateR` recommends, simply filter for
`recommended == TRUE`.

``` r
tp53_recommended <- filter(tp53_res$pairs, tp53_res$pairs$recommended == TRUE)
kableExtra::scroll_box(knitr::kable(head(tp53_recommended,5)),width = '700px')
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:700px; ">

<table>

<thead>

<tr>

<th style="text-align:right;">

upstream_pair
</th>

<th style="text-align:right;">

downstream_pair
</th>

<th style="text-align:right;">

exon_5p
</th>

<th style="text-align:right;">

exon_3p
</th>

<th style="text-align:left;">

compatible
</th>

<th style="text-align:left;">

frameshift
</th>

<th style="text-align:left;">

ptc_flag
</th>

<th style="text-align:left;">

terminal_exon_case
</th>

<th style="text-align:left;">

protospacer_sequence_5p
</th>

<th style="text-align:left;">

pam_sequence_5p
</th>

<th style="text-align:right;">

ontarget_score_5p
</th>

<th style="text-align:left;">

protospacer_sequence_3p
</th>

<th style="text-align:left;">

pam_sequence_3p
</th>

<th style="text-align:right;">

ontarget_score_3p
</th>

<th style="text-align:left;">

domains
</th>

<th style="text-align:left;">

recommended
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.7988059
</td>

<td style="text-align:left;">

ACTGGGACGGAACAGCTTTG
</td>

<td style="text-align:left;">

AGG
</td>

<td style="text-align:right;">

0.5867680
</td>

<td style="text-align:left;">

p53-like_TF_DNA-bd_sf; p53_tumour_suppressor; p53_DNA-bd;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.7988059
</td>

<td style="text-align:left;">

AGAGAATCTCCGCAAGAAAG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5479648
</td>

<td style="text-align:left;">

p53-like_TF_DNA-bd_sf; p53_tumour_suppressor; p53_DNA-bd;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

CCATTGTTCAATATCGTCCG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5255072
</td>

<td style="text-align:left;">

AGTGGTTTCTTCTTTGGCTG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5661123
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf; p53_tetrameristn; p53_tetramer_sf
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

CCTTCCCAGAAAACCTACCA
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5585539
</td>

<td style="text-align:left;">

AGTGGTTTCTTCTTTGGCTG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5661123
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf; p53_tetrameristn; p53_tetramer_sf
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

10
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

FALSE
</td>

<td style="text-align:left;">

CCCCGGACGATATTGAACAA
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

0.6092504
</td>

<td style="text-align:left;">

AGTGGTTTCTTCTTTGGCTG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

0.5661123
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf; p53_tetrameristn; p53_tetramer_sf
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

</tbody>

</table>

</div>

### Interpreting the `run_mutateR()` output

##### The `pairs` data.frame:

`upstream_exon` and `downstream_exon` refer to the two phase-compatible
exons, while `exon_5p` and `exon_3p` refer to the 5-prime and 3-prime
exons to simultaneously target using the selected gRNAs in that row of
the data.frame.

The `compatible` flag (`run_mutateR()` currently only returns compatible
exon pairs) highlights if the exon pair is phase-compatible; if a
non-compatible exon pair were used, one might expect the `frameshift`
and `ptc_flag` values to be `TRUE`. The `terminal_exon_case` flag is
used for a special case, where a frameshift has introduced a PTC into
the terminal exon of the selected transcript. Final-exon PTCs induce NMD
with far lower efficiency than those in non-terminal exons, and are
significantly less likely to result in transcriptional compensatory
behaviour.

The actual gRNA protospacer sequences are given in
`protospacer_sequence_5p` and `protospacer_sequence_3p` (with 5p and 3p
referring to the relative position in the transcript), and the
corresponding PAMs in `pam_sequence_5p` and `pam_sequence_3p`
respectively.

> Note that while both + and - strands are scanned for PAMs, the -
> strand sequences are automatically reverse-complemented so will always
> display NGG PAMs (instead of CCN).

Each 5-prime and 3-prime gRNA sequence also has corresponding
`ontarget_score`, `start` and `end` metadata. The `ontarget_score`
relates to whatever scoring method you selected when executing
`run_mutateR()`, while `start` and `end` are the genomic coordinates of
the start and end of the protospacer sequence.

`domains` contains protein domain annotations corresponding to the
selected exon pairs derived from the `map_protein_domains()` function,
while the `recommended` flag is `TRUE` if the gRNA pair meets all
validity criteria (including both gRNAs passing the on-target scoring
threshold).

##### Visualisation

`mutateR` also produces a basic visualisation for exon
phase-compatibility, by default a heatmap.

``` r
tp53_res$plot
```

![](README_files/figure-gfm/mutateR_plot_heatmap-1.png)<!-- -->

The x-axis represents the 5-prime exon, and the y-axis the three-prime
exon, to represent a given exon pair. UTR-containing exons are in red,
contiguous exons in green, while non-contiguous phase-compatible exons
are in yellow. All incompatible exon pairs are blue.

Domain annotations from Pfam are provided in a secondary plot element
below the heatmap x-axis to assist you in selecting regions should you
wish to ablate specific functional domains.

> A future release will allow interactive plotting to represent
> visualisation and easy retrieval of desired gRNA pairs directly from
> the plot.

The partially deprecated arc-based plotting mode can be accessed by
specifying `"arc"` in the `plot_mode` parameter in `run_mutateR()`. It
will remain, for now, the default plotting mode for small-gene edge
cases (described below).

This visualisation mode might be refined in future (if I feel like it),
but presently I will be focusing my effort on improving the heatmap.

``` r
tp53_arc <- run_mutateR(
  gene_id = "TP53",
  species = "hsapiens",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  nuclease = "Cas9",
  score_method = 'ruleset1',
  top_n = NULL,
  quiet = TRUE,
  plot_mode = 'arc'
  )
#> 
#> AGG CGG GGG TGG 
#> 102  36 136 125
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Flattening on-target scores from GRanges ...
#> Retrieving InterPro domain annotations from Ensembl Genes mart...
#> Generated 2664 candidate exon‑flanking gRNA pairs.
#> Plotting exon phase compatibility and gRNA pairs...

tp53_arc$plot
```

![](README_files/figure-gfm/mutateR_plot_arc-1.png)<!-- -->

##### Other data levels

The `run_mutateR()` output also contains two GRanges objects: `exons`
and `scored_grnas`.

`exons` contains the exon metadata (coordinates, strand, phase, length,
etc.):

``` r
head(tp53_res$exons,5)
#> GRanges object with 5 ranges and 8 metadata columns:
#>       seqnames          ranges strand | ensembl_exon_id      rank start_phase
#>          <Rle>       <IRanges>  <Rle> |     <character> <integer>   <integer>
#>   [1]       17 7687377-7687490      - | ENSE00003753508         1          -1
#>   [2]       17 7676521-7676622      - | ENSE00004023728         2          -1
#>   [3]       17 7676382-7676403      - | ENSE00002419584         3           2
#>   [4]       17 7675994-7676272      - | ENSE00003625790         4           0
#>   [5]       17 7675053-7675236      - | ENSE00003518480         5           0
#>       end_phase cds_start   cds_end transcript_cds_length exon_cds_length
#>       <integer> <integer> <integer>             <integer>       <integer>
#>   [1]        -1      <NA>      <NA>                  1182            <NA>
#>   [2]         2         1        74                  1182              74
#>   [3]         0        75        96                  1182              22
#>   [4]         0        97       375                  1182             279
#>   [5]         1       376       559                  1182             184
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

While `scored_grnas` contains the complete set of (unpaired) gRNAs,
their genomic coordinates, and the on-target score. The
`sequence_context` column is the protospacer + PAM as well as the
flanking nucleotides required for scoring (according to
`crisprScore::scoringMethodsInfo`) based on the selected scoring method.

``` r
head(tp53_res$scored_grnas,5)
#> GRanges object with 5 ranges and 8 metadata columns:
#>       seqnames          ranges strand | exon_rank protospacer_sequence
#>          <Rle>       <IRanges>  <Rle> | <integer>          <character>
#>   [1]    chr17 7687381-7687403      + |         1 AAAGTCTAGAGCCACCGTCC
#>   [2]    chr17 7687382-7687404      + |         1 AAGTCTAGAGCCACCGTCCA
#>   [3]    chr17 7687388-7687410      + |         1 AGAGCCACCGTCCAGGGAGC
#>   [4]    chr17 7687398-7687420      + |         1 TCCAGGGAGCAGGTAGCTGC
#>   [5]    chr17 7687399-7687421      + |         1 CCAGGGAGCAGGTAGCTGCT
#>       pam_sequence        target_sequence  cut_site       sequence_context
#>        <character>            <character> <numeric>            <character>
#>   [1]          AGG AAAGTCTAGAGCCACCGTCC..   7687398 CTCAAAAGTCTAGAGCCACC..
#>   [2]          GGG AAGTCTAGAGCCACCGTCCA..   7687399 TCAAAAGTCTAGAGCCACCG..
#>   [3]          AGG AGAGCCACCGTCCAGGGAGC..   7687405 GTCTAGAGCCACCGTCCAGG..
#>   [4]          TGG TCCAGGGAGCAGGTAGCTGC..   7687415 ACCGTCCAGGGAGCAGGTAG..
#>   [5]          GGG CCAGGGAGCAGGTAGCTGCT..   7687416 CCGTCCAGGGAGCAGGTAGC..
#>              gc ontarget_score
#>       <numeric>      <numeric>
#>   [1]  0.566667      0.0218384
#>   [2]  0.566667      0.1859041
#>   [3]  0.633333      0.1229175
#>   [4]  0.666667      0.0214385
#>   [5]  0.700000      0.2763628
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### Special cases

#### Small genes

Given that NMD induction requires persistence of EJCs, genes with one or
two exons tend to not be quite as prone to this mode of transcriptional
adaptive response. However, we still want these types of genes to be
able to be processed by `mutateR`, which will automatically recognise
these edge cases if you provide such a gene (or a specific transcript
with only 1-2 exons).

Instead of computing exon phase compatibility (which isn’t a factor for
these types of gene), the pipeline instead prioritises providing the
user with gRNA pairs that 1) score well per the applied on-target
method, and 2) produce a known (approximate) deletion size.

Here’s an example with zebrafish *opn4.1*, which contains a single,
large coding exon):

``` r
library(BSgenome.Drerio.UCSC.danRer11)

opn4.1 <- run_mutateR(gene_id = 'opn4.1',
                     species = 'drerio',
                     genome = BSgenome.Drerio.UCSC.danRer11,
                     nuclease = 'Cas9',
                     score_method = 'ruleset1')
#> Retrieving gene/transcript information...
#> Using transcript: ENSDART00000018501 for gene: opn4.1
#> Locating Cas9 target sites...
#> Scoring gRNAs using model: ruleset1
#> Computing on‑target scores using ruleset1 model...
#> Re‑orienting 129 guides so PAM = NGG at positions 26–27.
#> PAM triplet distribution (positions 25–27):
#> 
#> AGG CGG GGG TGG 
#>  71  43  68 105
#> Scored 287 guides using ruleset1.
#> Assembling valid gRNA pairs for opn4.1 ...
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Single-exon/two-exon gene detected: constructing intragenic deletion pairs.
#> Generated 41041 intragenic deletion pairs; 435 meet score cutoff.
#> Returning 435 recommended intragenic pairs.
#> Detected intragenic assembly mode (≤2 exons).
#> Plotting exon phase compatibility and gRNA pairs...
#> mutateR pipeline completed for opn4.1, finding 435 gRNA pairs.
```

``` r
opn4.1_recommended <- filter(opn4.1$pairs, opn4.1$pairs$recommended == TRUE)
kableExtra::scroll_box(knitr::kable(head(opn4.1_recommended,5)),width = '700px')
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:700px; ">

<table>

<thead>

<tr>

<th style="text-align:left;">

protospacer_sequence_5p
</th>

<th style="text-align:left;">

protospacer_sequence_3p
</th>

<th style="text-align:right;">

ontarget_score_5p
</th>

<th style="text-align:right;">

ontarget_score_3p
</th>

<th style="text-align:right;">

exon_5p
</th>

<th style="text-align:right;">

exon_3p
</th>

<th style="text-align:right;">

del_size
</th>

<th style="text-align:left;">

transcript_id
</th>

<th style="text-align:left;">

recommended
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

TTAATCAGGGGTCACCACAG
</td>

<td style="text-align:left;">

GGTGCAGTTGATGTCTCCAG
</td>

<td style="text-align:right;">

0.5483174
</td>

<td style="text-align:right;">

0.5378872
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

2777
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

ACATATGCTGGATAATTTGG
</td>

<td style="text-align:left;">

GGTGCAGTTGATGTCTCCAG
</td>

<td style="text-align:right;">

0.6630494
</td>

<td style="text-align:right;">

0.5378872
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

2768
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

CCAGCATATGTTTTACGCAA
</td>

<td style="text-align:left;">

GGTGCAGTTGATGTCTCCAG
</td>

<td style="text-align:right;">

0.6543132
</td>

<td style="text-align:right;">

0.5378872
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

2736
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

ACACACATATACACTACAGG
</td>

<td style="text-align:left;">

GGTGCAGTTGATGTCTCCAG
</td>

<td style="text-align:right;">

0.6773641
</td>

<td style="text-align:right;">

0.5378872
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

2661
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

TTAATCAGGGGTCACCACAG
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.5483174
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

2645
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

</tr>

</tbody>

</table>

</div>

There is a basic plotting mode currently implemented for this type of
gene, though this is mainly a placeholder and doesn’t really convey any
useful information yet.

#### Terminal exon PTCs

If a PTC occurs in the terminal exon of a gene, the PTC-bearing
transcript is unlikely to induce NMD and subsequent compensation (an
example of such an allele is the *sgsh<sup>Δex5−6</sup>* zebrafish
mutant in <a href="https://www.mdpi.com/1422-0067/22/11/5948"
style="font-size: 11pt;"
title="An Engineered sgsh Mutant Zebrafish Recapitulates Molecular and Behavioural Pathobiology of Sanfilippo Syndrome A/MPS IIIA">Douek
et al., IJMS (2021)</a>). `mutateR` factors in predicted terminal exon
PTCs when scanning for tolerated deletions, and flags these cases in the
`run_mutateR` output in the `$pairs` dataframe, under the
`terminal_exon_case` column.

## Manual function execution

You may wish to call lower-level functions directly (e.g., for use in
other programmatic pipelines). Note, however, that many of these
functions are built to take outputs from other `mutateR` functions.

### `get_gene_info()`

``` r
tp53_geneinfo <- get_gene_info(gene_id = 'TP53',
                               species = 'hsapiens',
                               id_type = 'symbol')
tp53_geneinfo
#> $canonical
#>    ensembl_gene_id ensembl_transcript_id chromosome_name start_position
#> 25 ENSG00000141510       ENST00000269305              17        7661779
#>    end_position strand external_gene_name transcript_is_canonical
#> 25      7687546      -               TP53                    TRUE
#> 
#> $all
#>    ensembl_gene_id ensembl_transcript_id chromosome_name start_position
#> 25 ENSG00000141510       ENST00000269305              17        7661779
#> 4  ENSG00000141510       ENST00000359597              17        7661779
#> 1  ENSG00000141510       ENST00000413465              17        7661779
#> 13 ENSG00000141510       ENST00000420246              17        7661779
#> 21 ENSG00000141510       ENST00000445888              17        7661779
#> 14 ENSG00000141510       ENST00000455263              17        7661779
#> 19 ENSG00000141510       ENST00000503591              17        7661779
#> 5  ENSG00000141510       ENST00000504290              17        7661779
#> 6  ENSG00000141510       ENST00000504937              17        7661779
#> 32 ENSG00000141510       ENST00000505014              17        7661779
#> 18 ENSG00000141510       ENST00000508793              17        7661779
#> 22 ENSG00000141510       ENST00000509690              17        7661779
#> 7  ENSG00000141510       ENST00000510385              17        7661779
#> 20 ENSG00000141510       ENST00000514944              17        7661779
#> 33 ENSG00000141510       ENST00000571370              17        7661779
#> 31 ENSG00000141510       ENST00000574684              17        7661779
#> 29 ENSG00000141510       ENST00000576024              17        7661779
#> 23 ENSG00000141510       ENST00000604348              17        7661779
#> 11 ENSG00000141510       ENST00000610292              17        7661779
#> 15 ENSG00000141510       ENST00000610538              17        7661779
#> 8  ENSG00000141510       ENST00000610623              17        7661779
#> 9  ENSG00000141510       ENST00000618944              17        7661779
#> 10 ENSG00000141510       ENST00000619186              17        7661779
#> 24 ENSG00000141510       ENST00000619485              17        7661779
#> 12 ENSG00000141510       ENST00000620739              17        7661779
#> 16 ENSG00000141510       ENST00000622645              17        7661779
#> 2  ENSG00000141510       ENST00000635293              17        7661779
#> 3  ENSG00000141510       ENST00000714356              17        7661779
#> 17 ENSG00000141510       ENST00000714357              17        7661779
#> 26 ENSG00000141510       ENST00000714358              17        7661779
#> 30 ENSG00000141510       ENST00000714359              17        7661779
#> 27 ENSG00000141510       ENST00000714408              17        7661779
#> 28 ENSG00000141510       ENST00000714409              17        7661779
#> 34 ENSG00000141510       ENST00000905353              17        7661779
#> 35 ENSG00000141510       ENST00000923566              17        7661779
#> 36 ENSG00000141510       ENST00000923567              17        7661779
#> 37 ENSG00000141510       ENST00000923568              17        7661779
#> 38 ENSG00000141510       ENST00000923569              17        7661779
#> 39 ENSG00000141510       ENST00000949117              17        7661779
#>    end_position strand external_gene_name transcript_is_canonical
#> 25      7687546      -               TP53                    TRUE
#> 4       7687546      -               TP53                      NA
#> 1       7687546      -               TP53                      NA
#> 13      7687546      -               TP53                      NA
#> 21      7687546      -               TP53                      NA
#> 14      7687546      -               TP53                      NA
#> 19      7687546      -               TP53                      NA
#> 5       7687546      -               TP53                      NA
#> 6       7687546      -               TP53                      NA
#> 32      7687546      -               TP53                      NA
#> 18      7687546      -               TP53                      NA
#> 22      7687546      -               TP53                      NA
#> 7       7687546      -               TP53                      NA
#> 20      7687546      -               TP53                      NA
#> 33      7687546      -               TP53                      NA
#> 31      7687546      -               TP53                      NA
#> 29      7687546      -               TP53                      NA
#> 23      7687546      -               TP53                      NA
#> 11      7687546      -               TP53                      NA
#> 15      7687546      -               TP53                      NA
#> 8       7687546      -               TP53                      NA
#> 9       7687546      -               TP53                      NA
#> 10      7687546      -               TP53                      NA
#> 24      7687546      -               TP53                      NA
#> 12      7687546      -               TP53                      NA
#> 16      7687546      -               TP53                      NA
#> 2       7687546      -               TP53                      NA
#> 3       7687546      -               TP53                      NA
#> 17      7687546      -               TP53                      NA
#> 26      7687546      -               TP53                      NA
#> 30      7687546      -               TP53                      NA
#> 27      7687546      -               TP53                      NA
#> 28      7687546      -               TP53                      NA
#> 34      7687546      -               TP53                      NA
#> 35      7687546      -               TP53                      NA
#> 36      7687546      -               TP53                      NA
#> 37      7687546      -               TP53                      NA
#> 38      7687546      -               TP53                      NA
#> 39      7687546      -               TP53                      NA
```

### `get_exon_structures()`

``` r
tp53_canonical_tx <- tp53_geneinfo$canonical$ensembl_transcript_id[1] 
tp53_exons <- get_exon_structures(transcript_id = tp53_canonical_tx,
                                  species = 'hsapiens',
                                  output = 'GRanges')
tp53_exons
#> GRanges object with 11 ranges and 8 metadata columns:
#>        seqnames          ranges strand | ensembl_exon_id      rank start_phase
#>           <Rle>       <IRanges>  <Rle> |     <character> <integer>   <integer>
#>    [1]       17 7687377-7687490      - | ENSE00003753508         1          -1
#>    [2]       17 7676521-7676622      - | ENSE00004023728         2          -1
#>    [3]       17 7676382-7676403      - | ENSE00002419584         3           2
#>    [4]       17 7675994-7676272      - | ENSE00003625790         4           0
#>    [5]       17 7675053-7675236      - | ENSE00003518480         5           0
#>    [6]       17 7674859-7674971      - | ENSE00003723991         6           1
#>    [7]       17 7674181-7674290      - | ENSE00003712342         7           0
#>    [8]       17 7673701-7673837      - | ENSE00003725258         8           2
#>    [9]       17 7673535-7673608      - | ENSE00003786593         9           1
#>   [10]       17 7670609-7670715      - | ENSE00003545950        10           0
#>   [11]       17 7668421-7669690      - | ENSE00004023724        11           2
#>        end_phase cds_start   cds_end transcript_cds_length exon_cds_length
#>        <integer> <integer> <integer>             <integer>       <integer>
#>    [1]        -1      <NA>      <NA>                  1182            <NA>
#>    [2]         2         1        74                  1182              74
#>    [3]         0        75        96                  1182              22
#>    [4]         0        97       375                  1182             279
#>    [5]         1       376       559                  1182             184
#>    [6]         0       560       672                  1182             113
#>    [7]         2       673       782                  1182             110
#>    [8]         1       783       919                  1182             137
#>    [9]         0       920       993                  1182              74
#>   [10]         2       994      1100                  1182             107
#>   [11]        -1      1101      1182                  1182              82
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### `find_cas9_sites()`

``` r
tp53_cas9_sites <- find_cas9_sites(exon_gr = tp53_exons,
                                   genome = BSgenome.Hsapiens.UCSC.hg38,
                                   pam = 'NGG',
                                   protospacer_length = 20)
tp53_cas9_sites
#> GRanges object with 399 ranges and 6 metadata columns:
#>         seqnames          ranges strand | exon_rank protospacer_sequence
#>            <Rle>       <IRanges>  <Rle> | <integer>          <character>
#>     [1]    chr17 7687381-7687403      + |         1 AAAGTCTAGAGCCACCGTCC
#>     [2]    chr17 7687382-7687404      + |         1 AAGTCTAGAGCCACCGTCCA
#>     [3]    chr17 7687388-7687410      + |         1 AGAGCCACCGTCCAGGGAGC
#>     [4]    chr17 7687398-7687420      + |         1 TCCAGGGAGCAGGTAGCTGC
#>     [5]    chr17 7687399-7687421      + |         1 CCAGGGAGCAGGTAGCTGCT
#>     ...      ...             ...    ... .       ...                  ...
#>   [395]    chr17 7668434-7668456      - |        11 TGAACATGAGTTTTTTATGG
#>   [396]    chr17 7668433-7668455      - |        11 GAACATGAGTTTTTTATGGC
#>   [397]    chr17 7668430-7668452      - |        11 CATGAGTTTTTTATGGCGGG
#>   [398]    chr17 7668412-7668434      - |        11 GGAGGTAGACTGACCCTTTT
#>   [399]    chr17 7668404-7668426      - |        11 ACTGACCCTTTTTGGACTTC
#>         pam_sequence        target_sequence  cut_site       sequence_context
#>          <character>            <character> <numeric>            <character>
#>     [1]          AGG AAAGTCTAGAGCCACCGTCC..   7687398 CTCAAAAGTCTAGAGCCACC..
#>     [2]          GGG AAGTCTAGAGCCACCGTCCA..   7687399 TCAAAAGTCTAGAGCCACCG..
#>     [3]          AGG AGAGCCACCGTCCAGGGAGC..   7687405 GTCTAGAGCCACCGTCCAGG..
#>     [4]          TGG TCCAGGGAGCAGGTAGCTGC..   7687415 ACCGTCCAGGGAGCAGGTAG..
#>     [5]          GGG CCAGGGAGCAGGTAGCTGCT..   7687416 CCGTCCAGGGAGCAGGTAGC..
#>     ...          ...                    ...       ...                    ...
#>   [395]          CGG TGAACATGAGTTTTTTATGG..   7668451 GTCTTGAACATGAGTTTTTT..
#>   [396]          GGG GAACATGAGTTTTTTATGGC..   7668450 TCTTGAACATGAGTTTTTTA..
#>   [397]          AGG CATGAGTTTTTTATGGCGGG..   7668447 TGAACATGAGTTTTTTATGG..
#>   [398]          TGG GGAGGTAGACTGACCCTTTT..   7668429 GGCGGGAGGTAGACTGACCC..
#>   [399]          AGG ACTGACCCTTTTTGGACTTC..   7668421 GTAGACTGACCCTTTTTGGA..
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### `find_cas12a_sites()`

``` r
tp53_cas12a_sites <- find_cas12a_sites(exon_gr = tp53_exons,
                                       genome = BSgenome.Hsapiens.UCSC.hg38,
                                       pam = "TTTV",
                                       protospacer_length = 23)
tp53_cas12a_sites
#> GRanges object with 82 ranges and 6 metadata columns:
#>        seqnames          ranges strand | exon_rank   protospacer_sequence
#>           <Rle>       <IRanges>  <Rle> | <integer>            <character>
#>    [1]    chr17 7687434-7687460      + |         1 CGTTCGGGCTGGGAGCGTGC..
#>    [2]    chr17 7687458-7687484      + |         1 CACGACGGTGACACGCTTCC..
#>    [3]    chr17 7687435-7687461      - |         1 GCACGCTCCCAGCCCGAACG..
#>    [4]    chr17 7687411-7687437      - |         1 GTGTCCCCGGAGCCCAGCAG..
#>    [5]    chr17 7676581-7676607      - |         2 ATGTTTCCTGACTCAGAGGG..
#>    ...      ...             ...    ... .       ...                    ...
#>   [78]    chr17 7668613-7668639      - |        11 GCAAATGGAAGTCCTGGGTG..
#>   [79]    chr17 7668608-7668634      - |        11 TGGAAGTCCTGGGTGCTTCT..
#>   [80]    chr17 7668557-7668583      - |        11 GACCCAAAACCCAAAATGGC..
#>   [81]    chr17 7668549-7668575      - |        11 ACCCAAAATGGCAGGGGAGG..
#>   [82]    chr17 7668542-7668568      - |        11 ATGGCAGGGGAGGGAGAGAT..
#>        pam_sequence        target_sequence  cut_site       sequence_context
#>         <character>            <character> <numeric>            <character>
#>    [1]         TTTG TTTGCGTTCGGGCTGGGAGC..   7687456 ACACTTTGCGTTCGGGCTGG..
#>    [2]         TTTC TTTCCACGACGGTGACACGC..   7687480 GTGCTTTCCACGACGGTGAC..
#>    [3]         GAAA GAAAGCACGCTCCCAGCCCG..   7687439 GTGCTTTCCACGACGGTGAC..
#>    [4]         CAAA CAAAGTGTCCCCGGAGCCCA..   7687415 ACACTTTGCGTTCGGGCTGG..
#>    [5]         GAAA GAAAATGTTTCCTGACTCAG..   7676585 ACATTTTCAGACCTATGGAA..
#>    ...          ...                    ...       ...                    ...
#>   [78]         CAAA CAAAGCAAATGGAAGTCCTG..   7668617 TTGCTTTGTCCCGGGGCTCC..
#>   [79]         CAAA CAAATGGAAGTCCTGGGTGC..   7668612 TCCATTTGCTTTGTCCCGGG..
#>   [80]         CAAA CAAAGACCCAAAACCCAAAA..   7668561 GGTCTTTGAACCCTTGCTTG..
#>   [81]         CAAA CAAAACCCAAAATGGCAGGG..   7668553 GGGTTTTGGGTCTTTGAACC..
#>   [82]         CAAA CAAAATGGCAGGGGAGGGAG..   7668546 CCATTTTGGGTTTTGGGTCT..
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### `score_grnas()`

N.B. Currently only working for Cas9 and `ruleset1`. More rulesets and
nucleases to be implemented soon! Off-target analysis via `crisprScore`
for Cas9 (`mit`, `cfd`) also to be implemented in future.

``` r
tp53_cas9_gRNAs_scored <- score_grnas(grna_gr = tp53_cas9_sites,
                                      method = 'ruleset1')
#> Computing on‑target scores using ruleset1 model...
#> Re‑orienting 214 guides so PAM = NGG at positions 26–27.
#> PAM triplet distribution (positions 25–27):
#> 
#> AGG CGG GGG TGG 
#> 102  36 136 125
#> Scored 399 guides using ruleset1.
tp53_cas9_gRNAs_scored
#> GRanges object with 399 ranges and 8 metadata columns:
#>         seqnames          ranges strand | exon_rank protospacer_sequence
#>            <Rle>       <IRanges>  <Rle> | <integer>          <character>
#>     [1]    chr17 7687381-7687403      + |         1 AAAGTCTAGAGCCACCGTCC
#>     [2]    chr17 7687382-7687404      + |         1 AAGTCTAGAGCCACCGTCCA
#>     [3]    chr17 7687388-7687410      + |         1 AGAGCCACCGTCCAGGGAGC
#>     [4]    chr17 7687398-7687420      + |         1 TCCAGGGAGCAGGTAGCTGC
#>     [5]    chr17 7687399-7687421      + |         1 CCAGGGAGCAGGTAGCTGCT
#>     ...      ...             ...    ... .       ...                  ...
#>   [395]    chr17 7668434-7668456      - |        11 TGAACATGAGTTTTTTATGG
#>   [396]    chr17 7668433-7668455      - |        11 GAACATGAGTTTTTTATGGC
#>   [397]    chr17 7668430-7668452      - |        11 CATGAGTTTTTTATGGCGGG
#>   [398]    chr17 7668412-7668434      - |        11 GGAGGTAGACTGACCCTTTT
#>   [399]    chr17 7668404-7668426      - |        11 ACTGACCCTTTTTGGACTTC
#>         pam_sequence        target_sequence  cut_site       sequence_context
#>          <character>            <character> <numeric>            <character>
#>     [1]          AGG AAAGTCTAGAGCCACCGTCC..   7687398 CTCAAAAGTCTAGAGCCACC..
#>     [2]          GGG AAGTCTAGAGCCACCGTCCA..   7687399 TCAAAAGTCTAGAGCCACCG..
#>     [3]          AGG AGAGCCACCGTCCAGGGAGC..   7687405 GTCTAGAGCCACCGTCCAGG..
#>     [4]          TGG TCCAGGGAGCAGGTAGCTGC..   7687415 ACCGTCCAGGGAGCAGGTAG..
#>     [5]          GGG CCAGGGAGCAGGTAGCTGCT..   7687416 CCGTCCAGGGAGCAGGTAGC..
#>     ...          ...                    ...       ...                    ...
#>   [395]          CGG TGAACATGAGTTTTTTATGG..   7668451 GTCTTGAACATGAGTTTTTT..
#>   [396]          GGG GAACATGAGTTTTTTATGGC..   7668450 TCTTGAACATGAGTTTTTTA..
#>   [397]          AGG CATGAGTTTTTTATGGCGGG..   7668447 TGAACATGAGTTTTTTATGG..
#>   [398]          TGG GGAGGTAGACTGACCCTTTT..   7668429 GGCGGGAGGTAGACTGACCC..
#>   [399]          AGG ACTGACCCTTTTTGGACTTC..   7668421 GTAGACTGACCCTTTTTGGA..
#>                gc ontarget_score
#>         <numeric>      <numeric>
#>     [1]  0.566667      0.0218384
#>     [2]  0.566667      0.1859041
#>     [3]  0.633333      0.1229175
#>     [4]  0.666667      0.0214385
#>     [5]  0.700000      0.2763628
#>     ...       ...            ...
#>   [395]  0.433333      0.4110943
#>   [396]  0.433333      0.0695857
#>   [397]  0.433333      0.2398178
#>   [398]  0.566667      0.0744060
#>   [399]  0.500000      0.0161947
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## To be implemented

- More nucleases (maybe from CasPEDIA or similar?)

- All `crisprScore` scoring methods

- Visualisation improvements – in progress

- Workflow for handling of intronic target sequences

## Session information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26100)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_Australia.utf8  LC_CTYPE=English_Australia.utf8   
#> [3] LC_MONETARY=English_Australia.utf8 LC_NUMERIC=C                      
#> [5] LC_TIME=English_Australia.utf8    
#> 
#> time zone: Australia/Sydney
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] BSgenome.Drerio.UCSC.danRer11_1.4.2 RColorBrewer_1.1-3                 
#>  [3] tibble_3.3.0                        ggplot2_4.0.0                      
#>  [5] purrr_1.1.0                         dplyr_1.1.4                        
#>  [7] BSgenome.Hsapiens.UCSC.hg38_1.4.5   BSgenome_1.77.2                    
#>  [9] rtracklayer_1.69.1                  BiocIO_1.19.0                      
#> [11] Biostrings_2.77.2                   XVector_0.49.1                     
#> [13] GenomicRanges_1.61.5                GenomeInfoDb_1.45.12               
#> [15] Seqinfo_0.99.2                      IRanges_2.43.2                     
#> [17] S4Vectors_0.47.2                    BiocGenerics_0.55.1                
#> [19] generics_0.1.4                      mutateR_0.0.1                      
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.2.3                   bitops_1.0-9               
#>  [3] httr2_1.2.1                 biomaRt_2.65.16            
#>  [5] rlang_1.1.6                 magrittr_2.0.4             
#>  [7] matrixStats_1.5.0           compiler_4.5.1             
#>  [9] RSQLite_2.4.3               dir.expiry_1.17.0          
#> [11] png_0.1-8                   systemfonts_1.3.1          
#> [13] vctrs_0.6.5                 stringr_1.5.2              
#> [15] pkgconfig_2.0.3             crayon_1.5.3               
#> [17] fastmap_1.2.0               dbplyr_2.5.1               
#> [19] labeling_0.4.3              Rsamtools_2.25.3           
#> [21] rmarkdown_2.30              UCSC.utils_1.5.0           
#> [23] bit_4.6.0                   xfun_0.54                  
#> [25] randomForest_4.7-1.2        cachem_1.1.0               
#> [27] jsonlite_2.0.0              progress_1.2.3             
#> [29] blob_1.2.4                  DelayedArray_0.35.3        
#> [31] BiocParallel_1.43.4         parallel_4.5.1             
#> [33] prettyunits_1.2.0           R6_2.6.1                   
#> [35] stringi_1.8.7               reticulate_1.44.0          
#> [37] Rcpp_1.1.0                  SummarizedExperiment_1.39.2
#> [39] knitr_1.50                  Matrix_1.7-4               
#> [41] tidyselect_1.2.1            rstudioapi_0.17.1          
#> [43] dichromat_2.0-0.1           abind_1.4-8                
#> [45] yaml_2.3.10                 codetools_0.2-20           
#> [47] curl_7.0.0                  lattice_0.22-7             
#> [49] Biobase_2.69.1              withr_3.0.2                
#> [51] KEGGREST_1.49.1             S7_0.2.0                   
#> [53] evaluate_1.0.5              crisprScoreData_1.13.0     
#> [55] BiocFileCache_2.99.6        xml2_1.4.1                 
#> [57] ExperimentHub_2.99.5        pillar_1.11.1              
#> [59] BiocManager_1.30.26         filelock_1.0.3             
#> [61] MatrixGenerics_1.21.0       crisprScore_1.13.1         
#> [63] RCurl_1.98-1.17             BiocVersion_3.22.0         
#> [65] hms_1.1.4                   scales_1.4.0               
#> [67] glue_1.8.0                  tools_4.5.1                
#> [69] AnnotationHub_3.99.6        GenomicAlignments_1.45.4   
#> [71] XML_3.99-0.19               grid_4.5.1                 
#> [73] AnnotationDbi_1.71.1        basilisk_1.21.5            
#> [75] restfulr_0.0.16             cli_3.6.5                  
#> [77] rappdirs_0.3.3              kableExtra_1.4.0           
#> [79] textshaping_1.0.4           viridisLite_0.4.2          
#> [81] S4Arrays_1.9.1              svglite_2.2.2              
#> [83] gtable_0.3.6                digest_0.6.37              
#> [85] SparseArray_1.9.1           rjson_0.2.23               
#> [87] farver_2.1.2                memoise_2.0.1              
#> [89] htmltools_0.5.8.1           lifecycle_1.0.4            
#> [91] httr_1.4.7                  bit64_4.6.0-1
```

## References

El-Brolosy, M. A. *et al.* Genetic compensation triggered by mutant mRNA
degradation. *Nature* **568**, 193–197 (2019).
<https://doi.org/10.1038/s41586-019-1064-z>

Hoberecht, L., Perampalam, P., Lun, A. & Fortin, J. P. A comprehensive
Bioconductor ecosystem for the design of CRISPR guide RNAs across
nucleases and technologies. *Nat Commun* **13**, 6568 (2022).
<https://doi.org/10.1038/s41467-022-34320-7>

Lindeboom, R. G., Supek, F. & Lehner, B. The rules and impact of
nonsense-mediated mRNA decay in human cancers. *Nat Genet* **48**,
1112–1118 (2016). <https://doi.org/10.1038/ng.3664>

Ma, Z. *et al.* PTC-bearing mRNA elicits a genetic compensation response
via Upf3a and COMPASS components. *Nature* **568**, 259–263 (2019).
<https://doi.org/10.1038/s41586-019-1057-y>

Douek, A. M. *et al.* An Engineered sgsh Mutant Zebrafish Recapitulates
Molecular and Behavioural Pathobiology of Sanfilippo Syndrome A/MPS
IIIA. *Int J Mol Sci* **22**, 5948 (2021).
<https://doi.org/10.3390/ijms22115948>

Scott, D. A. & Zhang, F. Implications of human genetic variation in
CRISPR-based therapeutic genome editing. *Nat Med* **23**, 1095–1101
(2017). <https://doi.org/10.1038/nm.4377>

Yang, W., Zhu, J. K. & Jin, W. A catalog of gene editing sites and
genetic variations in editing sites in model organisms. *BMC Genomics*
**25**, 1153 (2024). <https://doi.org/10.1186/s12864-024-11073-9>
