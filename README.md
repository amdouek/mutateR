mutateR
================
Alon M Douek
2025-12-06

- [Overview](#overview)
  - [Background](#background)
  - [Why does this matter?](#why-does-this-matter)
  - [How can we reduce the likelihood of transcriptional
    adaptation?](#how-can-we-reduce-the-likelihood-of-transcriptional-adaptation)
- [How does `mutateR` work?](#how-does-mutater-work)
- [Getting started](#getting-started)
  - [Installation](#installation)
  - [Setting up the `mutateR` Python
    backend](#setting-up-the-mutater-python-backend)
- [The `mutateR` workflow](#the-mutater-workflow)
  - [Interpreting the `run_mutateR()`
    output](#interpreting-the-run_mutater-output)
  - [A note on primer design](#a-note-on-primer-design)
  - [Special cases](#special-cases)
- [Using the `mutateR_viewer`](#using-the-mutater_viewer)
- [Manual function execution](#manual-function-execution)
- [To be implemented](#to-be-implemented)
- [Session information](#session-information)
- [References](#references)

## Overview

`mutateR` streamlines design of frame-aware CRISPR/Cas-mediated
deletions.

It identifies pairs of non-contiguous exons in a given transcript that
are phase-compatible with each other (*i.e.*, pairs of exons that do not
normally neighbour each other in the transcript, but the 5’ exon ends in
the same reading frame position (0, 1, or 2) as the 3’ exon begins),
then detects gRNAs in each exon that facilitate their deletion such that
the resulting mutant allele preserves the original reading frame.

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

**An important caveat:** Transcriptional adaptation is not universal
(i.e., some genes are more prone to this behaviour than others, and it
is observed more frequently in some biological contexts than in others).
While `mutateR` assumes that your gene of interest is prone to
transcriptional adaptation when recommending gRNA pairs, you should
always base your mutagenesis rationale on the requirements of your
specific experimental design.

### How can we reduce the likelihood of transcriptional adaptation?

In short, **use tools for targeted mutagenesis to their strengths**.
Back when the only way to make a mutant animal was to randomly
mutagenise its genome and then recover alleles based on phenotypes,
researchers didn’t have a great deal of control over the nature of the
resulting alleles. Modern sequence-targeting methods (such as
CRISPR/Cas9) allow for much more control over the nature of the
mutations that we can generate.

If we delete parts of a gene that correspond to functional domains of
the encoded protein, but maintain the original reading frame (such that
there is no PTC), we are left with a stably-expressed, truncated mutant
transcript that does not encode a functional protein but also does not
induce transcriptional adaptation via NMD.

> *Note*: There is nuance to NMD induction. For example, PTCs are much
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
protospacers that do not cross exon boundaries.

After calculating on-target scores for potential gRNAs using rule sets
from `crisprScore` ([Hoberecht et al., *Nature Communications*
(2022)](https://www.nature.com/articles/s41467-022-34320-7 "A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies"))
where possible, or by manual re-implementation of methods not available
from `crisprScore`, `mutateR` returns recommended gRNA pairs to target
exons where the flanking exons are in-frame with each other.

It also produces a graphical representation of the selected transcript,
its phase-compatible exons, and the recommended exon pairs for
targeting. Successful `mutateR` pipeline executions can now be passed to
the `mutateR_viewer()` function for further inspection in an interactive
Shiny app.

> As of December 2025, `mutateR` utilises a Python backend (via
> `reticulate`) to allow for implementation of python-dependent
> functionalities. This was originally done to implement
> deep-learning-based on-target scoring methods (dependent on TensorFlow
> etc.) which could not be accessed by Windows OS users from
> `crisprScore`, but now is also required for the primer design facet of
> the `run_mutateR` pipeline (via `primer3-py`).
>
> As these core functionalities (and likely future ones) use the Python
> backend, it is **highly recommended** that all `mutateR` users install
> and activate the `r-mutater` Python env (containing all the requisite
> modules) from within `mutateR`. See the relevant Installation section
> (\[Setting up the mutateR Python backend\]) below for more
> information.

The `mutateR` package consists of the following ordered functions:

| Step | Function | Functionality |
|----|----|----|
| 1 | `get_gene_info()` | Query Ensembl for gene and transcript metadata. |
| 2 | `get_exon_structures()` | Retrieve exon coordinates and phase metadata. |
| 3 | `find_cas9_sites()`/`find_cas12a_sites()` | Identify nuclease recognition sites. |
| 4 | `check_exon_phase()` | Determine phase-compatible exon pairs. |
| 5 | `check_frameshift_ptc()` | Calculate and flag frameshift/PTC consequences. |
| 6 | `map_protein_domains()` | Retrieve protein domain annotations from InterPro corresponding to exons. |
| 7 | `score_grnas()` | Score guide RNAs via `crisprScore`. |
| 8 | `filter_valid_grnas()` and `assemble_grna_pairs()` | Filter allowed gRNA pairs and assemble into a readable data.frame. |
| 9 | `plot_grna_design()` | Visualise phase compatibility of non-contiguous exon pairs and top-ranked valid gRNA pairs. |
| 10 | `plot_grna_heatmap()` | Helper function for `plot_grna_design()` but can be executed separately. |
| 11 | `plot_grna_interactive()` | Helper function for interactive plotting mode. |
| 12 | `run_mutateR()` | Wrapper; runs entire pipeline. |
| 13 | `mutateR_viewer()` | Launches an Shiny app for gRNA extraction from an interactive heatmap. |

`mutateR` also has an additional function, `design_grna_pairs()`, that
is not included in the wrapper function. This function can be used in
isolation to produce a data.frame of gRNA pairs for the gene of
interest. Other scripts such as `deep_learning_utils.R`, `python_env.R`,
`primer3-backend.R` and `design_primers.R` contain helper functions that
are critical but not directly user-executed pipeline components.

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

### Setting up the `mutateR` Python backend

In order to use the full suite of functionalities in `mutateR` (e.g.,
automated genotyping primer design, all gRNA on-target scoring methods),
you must install and activate the package’s Conda environment. This env
contains all the necessary dependencies for Python-reliant package
functions.

#### Install `reticulate`

`reticulate` is the [R interface to
Python](https://rstudio.github.io/reticulate/) required to run `mutateR`
as intended. You can install it from CRAN if needed, but it is imported
upon installation of `mutateR`.

#### Install the `r-mutater` Python env

Installing the Python backend is simple: In R, after installing and
loading the `mutateR` package, simply run the command
`install_mutater_env()`.

> If you have never used Python before, it likely does not yet exist on
> your device. If this is the case, and you are running `mutateR` for
> the first time, this command will install
> [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
> via `reticulate`. This might take a few minutes - ensure it completes
> fully before proceeding.

This command will create an env called `r-mutater` with the specified
`python_version` containing the following dependencies:
`"tensorflow-cpu", "numpy<2", "h5py", "pandas", "scipy", "primer3-py"`.

If for whatever reason the env has not been installed correctly (or if a
new dependency has been added), rerun the installation with the
parameter `fresh = TRUE` to remove and reinstall the env.

**Important:** Restart your R session (Ctrl/Cmd + Shift + F10) after
installation is complete.

#### Activate the env

After **restarting your R session**, run `activate_mutater_env()` to
point your R session to the `r-mutater` env. `reticulate`’s default env
is `RETICULATE_PYTHON`; `activate_mutater_env()` will unset this and set
it as `r-mutater`. This will only last for the current R session and
will revert to the default upon session restart (i.e., you will need to
reactivate the env each time you initiate a new session, but you will
only need to *reinstall* the env if a new dependency is added or if it
is corrupted).

Once the env is successfully activated, `[1] TRUE` will be printed in
the console. At this point, you will be able to make use of the full
suite of `mutateR` functionalities.

## The `mutateR` workflow

In the below example run, we will run `mutateR` on human *TP53* to
select optimal Cas9 gRNA pairs*.* Assuming you have installed both
`mutateR` and your BSgenome of choice (in this case, human), and set up
the Python backend correctly:

``` r
library(reticulate)
library(mutateR)
library(BSgenome.Hsapiens.UCSC.hg38)

activate_mutater_env()
#> [1] TRUE
```

The wrapper function `run_mutateR()` executes the full pipeline (using
the canonical transcript by default). For human *TP53:*

``` r
tp53_res <- run_mutateR(
  gene_id = "TP53",
  species = "hsapiens",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  nuclease = "Cas9",
  score_method = 'ruleset1',
  design_primers = TRUE,
  quiet = FALSE,
  plot_mode = 'heat'
  )
#> Retrieving gene/transcript information...
#> Using transcript: ENST00000269305 for gene: TP53
#> Locating Cas9 target sites...
#> PAM distribution:
#> 
#> AGG CGG GGG TGG 
#> 104  39 133 123
#> Scoring gRNAs using model: ruleset1
#> Computing on‑target scores using ruleset1 model...
#> Scored 399 guides using ruleset1.
#> Assembling valid gRNA pairs for TP53 ...
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Detected probability-based scores (e.g. ruleset1). Using default cutoff: 0.5
#> Retrieving InterPro domain annotations from Ensembl Genes mart...
#> Generated 2784 candidate exon‑flanking gRNA pairs.
#> Designing genotyping primers for 35 recommended pairs...
#> Preparing primer design batch requests...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   9%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  26%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  34%  |                                                                              |==========================                                            |  37%  |                                                                              |============================                                          |  40%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  60%  |                                                                              |============================================                          |  63%  |                                                                              |==============================================                        |  66%  |                                                                              |================================================                      |  69%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  77%  |                                                                              |========================================================              |  80%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  91%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
#> Running Batch Primer3 (35 designs)...
#> Mapping results to dataframe...
#> Designed primers for 35/35 pairs.
#> Plotting exon phase compatibility and gRNA pairs...
#> Retrieving Pfam domain annotations from Ensembl Genes mart...
#> mutateR pipeline completed for TP53, finding 2784 gRNA pairs.
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

<th style="text-align:left;">

seqnames_5p
</th>

<th style="text-align:right;">

start_5p
</th>

<th style="text-align:right;">

end_5p
</th>

<th style="text-align:left;">

protospacer_sequence_5p
</th>

<th style="text-align:left;">

pam_sequence_5p
</th>

<th style="text-align:right;">

cut_site_5p
</th>

<th style="text-align:right;">

ontarget_score_5p
</th>

<th style="text-align:left;">

seqnames_3p
</th>

<th style="text-align:right;">

start_3p
</th>

<th style="text-align:right;">

end_3p
</th>

<th style="text-align:left;">

protospacer_sequence_3p
</th>

<th style="text-align:left;">

pam_sequence_3p
</th>

<th style="text-align:right;">

cut_site_3p
</th>

<th style="text-align:right;">

ontarget_score_3p
</th>

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

<th style="text-align:right;">

genomic_deletion_size
</th>

<th style="text-align:right;">

transcript_deletion_size
</th>

<th style="text-align:left;">

domains
</th>

<th style="text-align:left;">

recommended
</th>

<th style="text-align:left;">

priming_strategy
</th>

<th style="text-align:left;">

primer_ext_fwd
</th>

<th style="text-align:left;">

primer_ext_rev
</th>

<th style="text-align:left;">

primer_int_fwd
</th>

<th style="text-align:left;">

primer_int_rev
</th>

<th style="text-align:left;">

exp_wt_size
</th>

<th style="text-align:right;">

exp_mut_size
</th>

<th style="text-align:right;">

exp_int_size
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676034
</td>

<td style="text-align:right;">

7676056
</td>

<td style="text-align:left;">

CCCCGGACGATATTGAACAA
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7676051
</td>

<td style="text-align:right;">

0.6092504
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1165
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1424
</td>

<td style="text-align:right;">

259
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676080
</td>

<td style="text-align:right;">

7676102
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

AGG
</td>

<td style="text-align:right;">

7676097
</td>

<td style="text-align:right;">

0.6824355
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1211
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1424
</td>

<td style="text-align:right;">

213
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676189
</td>

<td style="text-align:right;">

7676211
</td>

<td style="text-align:left;">

CCTTCCCAGAAAACCTACCA
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

7676206
</td>

<td style="text-align:right;">

0.5585539
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1320
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

568
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676204
</td>

<td style="text-align:right;">

7676226
</td>

<td style="text-align:left;">

CAGAATGCAAGAAGCCCAGA
</td>

<td style="text-align:left;">

CGG
</td>

<td style="text-align:right;">

7676221
</td>

<td style="text-align:right;">

0.6499000
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1335
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

553
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676151
</td>

<td style="text-align:right;">

7676173
</td>

<td style="text-align:left;">

GAAGGGACAGAAGATGACAG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

7676168
</td>

<td style="text-align:right;">

0.7889533
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1282
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

606
</td>

<td style="text-align:right;">

NA
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

<th style="text-align:left;">

seqnames_5p
</th>

<th style="text-align:right;">

start_5p
</th>

<th style="text-align:right;">

end_5p
</th>

<th style="text-align:left;">

protospacer_sequence_5p
</th>

<th style="text-align:left;">

pam_sequence_5p
</th>

<th style="text-align:right;">

cut_site_5p
</th>

<th style="text-align:right;">

ontarget_score_5p
</th>

<th style="text-align:left;">

seqnames_3p
</th>

<th style="text-align:right;">

start_3p
</th>

<th style="text-align:right;">

end_3p
</th>

<th style="text-align:left;">

protospacer_sequence_3p
</th>

<th style="text-align:left;">

pam_sequence_3p
</th>

<th style="text-align:right;">

cut_site_3p
</th>

<th style="text-align:right;">

ontarget_score_3p
</th>

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

<th style="text-align:right;">

genomic_deletion_size
</th>

<th style="text-align:right;">

transcript_deletion_size
</th>

<th style="text-align:left;">

domains
</th>

<th style="text-align:left;">

recommended
</th>

<th style="text-align:left;">

priming_strategy
</th>

<th style="text-align:left;">

primer_ext_fwd
</th>

<th style="text-align:left;">

primer_ext_rev
</th>

<th style="text-align:left;">

primer_int_fwd
</th>

<th style="text-align:left;">

primer_int_rev
</th>

<th style="text-align:left;">

exp_wt_size
</th>

<th style="text-align:right;">

exp_mut_size
</th>

<th style="text-align:right;">

exp_int_size
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676034
</td>

<td style="text-align:right;">

7676056
</td>

<td style="text-align:left;">

CCCCGGACGATATTGAACAA
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7676051
</td>

<td style="text-align:right;">

0.6092504
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1165
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1424
</td>

<td style="text-align:right;">

259
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676080
</td>

<td style="text-align:right;">

7676102
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

AGG
</td>

<td style="text-align:right;">

7676097
</td>

<td style="text-align:right;">

0.6824355
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1211
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGAAGCTCCCAGAATGCCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1424
</td>

<td style="text-align:right;">

213
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676189
</td>

<td style="text-align:right;">

7676211
</td>

<td style="text-align:left;">

CCTTCCCAGAAAACCTACCA
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

7676206
</td>

<td style="text-align:right;">

0.5585539
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1320
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

568
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676204
</td>

<td style="text-align:right;">

7676226
</td>

<td style="text-align:left;">

CAGAATGCAAGAAGCCCAGA
</td>

<td style="text-align:left;">

CGG
</td>

<td style="text-align:right;">

7676221
</td>

<td style="text-align:right;">

0.6499000
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1335
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

553
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7676151
</td>

<td style="text-align:right;">

7676173
</td>

<td style="text-align:left;">

GAAGGGACAGAAGATGACAG
</td>

<td style="text-align:left;">

GGG
</td>

<td style="text-align:right;">

7676168
</td>

<td style="text-align:right;">

0.7889533
</td>

<td style="text-align:left;">

chr17
</td>

<td style="text-align:right;">

7674869
</td>

<td style="text-align:right;">

7674891
</td>

<td style="text-align:left;">

TCCTCAGCATCTTATCCGAG
</td>

<td style="text-align:left;">

TGG
</td>

<td style="text-align:right;">

7674886
</td>

<td style="text-align:right;">

0.7988059
</td>

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

<td style="text-align:right;">

1282
</td>

<td style="text-align:right;">

576
</td>

<td style="text-align:left;">

p53_tumour_suppressor; p53_DNA-bd; p53-like_TF_DNA-bd_sf; p53_TAD2;
p53/RUNT-type_TF_DNA-bd_sf
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

TAAGCAGCAGGAGAAAGCCC
</td>

<td style="text-align:left;">

TGCTGGATCCCCACTTTTCC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1888
</td>

<td style="text-align:right;">

606
</td>

<td style="text-align:right;">

NA
</td>

</tr>

</tbody>

</table>

</div>

### Interpreting the `run_mutateR()` output

#### The `pairs` data.frame

`$pairs` contains several columns - some of these are not strictly
useful for the end-user (e.g. genomic coordinates), but most represent
the primary output of the pipeline. The below sections will elaborate on
these:

##### Flanking and targeted exons

`upstream_exon` and `downstream_exon` refer to the two phase-compatible
exons, while `exon_5p` and `exon_3p` refer to the 5-prime and 3-prime
exons to simultaneously target using the selected gRNAs in that row of
the data.frame.

##### Pair classification flags

Each gRNA pair is annotated by several classifications. The `compatible`
flag (`run_mutateR()` currently only returns compatible exon pairs)
highlights if the exon pair is phase-compatible; if a non-compatible
exon pair were used, one might expect the `frameshift` and `ptc_flag`
values to be `TRUE`. The `terminal_exon_case` flag is used for a special
case, where a frameshift has introduced a PTC into the terminal exon of
the selected transcript. Final-exon PTCs induce NMD with far lower
efficiency than those in non-terminal exons, and are significantly less
likely to result in transcriptional compensatory behaviour.

##### gRNA pairs and associated deletion sizes

The actual gRNA protospacer sequences are given in
`protospacer_sequence_5p` and `protospacer_sequence_3p` (with 5p and 3p
referring to the relative position in the transcript of the gRNAs in the
given pair), and the corresponding PAMs in `pam_sequence_5p` and
`pam_sequence_3p` respectively.

> Note that while both + and - strands are scanned for PAMs, the -
> strand sequences are automatically reverse-complemented so will always
> display NGG PAMs (instead of CCN).

For each gRNA pair, the expected deletion size at both the genomic and
transcript level is reported in `genomic_deletion_size` and
`transcript_deletion_size`, respectively.

##### Genotyping primers and amplicon sizes

`run_mutateR()` now also automatically generates a pair (or two,
depending on the nature of the deletion - see [A note on primer
design](#a-note-on-primer-design) below) of primers for genotyping the
resulting deletion produced by the gRNA targeting. `priming_strategy`
reports the genotyping approach used (either ‘Flanking’ \[one primer
pair\] or ‘Dual Pair’ \[two primer pairs\]); `primer_ext_fwd` and
`primer_ext_rev` for the deletion-flanking primer sequences, and
`primer_int_fwd` and `primer_int_rev` for primer pair sequences nested
within the deleted region (Dual Pair only). `exp_wt_size` reports the
expected wild type amplicon size, and `exp_mut_size` the expected mutant
amplicon size from the flanking primers. `exp_int_size` reports the
amplicon size for the nested primer pair.

##### gRNA scoring

The currently working scoring methods are `"ruleset1"` and
`"deepspcas9"` for Cas9, and `"deepcpf1"` for Cas12a. More will be added
over time!

`ruleset1` is implemented via `crisprScore` and works in an R-only
environment, but the latter two require activation of the `r-mutater`
Python env. The weights for the deep learning-based methods can be found
in `inst/extdata`. Please note that in the context of `mutateR`,
`deepcpf1` refers to the
[Seq-DeepCpf1](https://www.nature.com/articles/nbt.4061) model
(sequence-only, without training on chromatin accessibility
information).

> Please also note that the authors of
> [DeepSpCas9](https://www.science.org/doi/10.1126/sciadv.aax9249)
> provided the training data (as TensorFlow checkpoint files) rather
> than as explicit weights. I have attempted to reconstruct the weights
> based on the model topology reported in the paper, and while it seems
> to be reporting values within the range expected, these scores should
> be treated with caution until I can do some more validation that my
> weights are equivalent to the ones used by the authors.

Each 5-prime and 3-prime gRNA sequence also has corresponding
`ontarget_score` relating to whatever scoring method you selected when
executing `run_mutateR()`. Note that the nature of the score will differ
between methods (e.g., `"ruleset1"` reports a probability value between
0-1, while `"deepspcas9"` reports a scalar linear regression value
(usually 0-100, occasionally negative).

##### Other

`domains` contains protein domain annotations corresponding to the
region to be deleted derived from the `map_protein_domains()` function,
while the `recommended` flag is `TRUE` if the gRNA pair meets all
validity criteria (including both gRNAs passing the on-target scoring
threshold). Importantly, genotyping primers are only generated for
recommended gRNA pairs.

#### Visualisation

##### Static visualisation modes

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
  design_primers = FALSE,
  top_n = NULL,
  quiet = TRUE,
  plot_mode = 'arc'
  )
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Detected probability-based scores (e.g. ruleset1). Using default cutoff: 0.5
#> Retrieving InterPro domain annotations from Ensembl Genes mart...
#> Generated 2784 candidate exon‑flanking gRNA pairs.
#> Plotting exon phase compatibility and gRNA pairs...

tp53_arc$plot
```

![](README_files/figure-gfm/mutateR_plot_arc-1.png)<!-- -->

##### Interactive visualisation modes

`mutateR` now offers an interactive heatmap plotting mode which can be
evoked by specifying `interactive = TRUE` inside `run_mutateR()`.
Hovering over a cell in the heatmap will report metadata for that exon
pair in a tooltip (e.g., whether the exon pair is (in)compatible, how
many gRNA pairs exist for that exon pair, etc.). This visualisation mode
is used as the basis for the much more useful `mutateR_viewer` Shiny app
(see below).

![](README_files/figure-gfm/heat_int.gif)

#### Other data levels

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
#> GRanges object with 5 ranges and 9 metadata columns:
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
#>              gc ontarget_score scoring_method
#>       <numeric>      <numeric>    <character>
#>   [1]  0.566667      0.0218384       ruleset1
#>   [2]  0.566667      0.1859041       ruleset1
#>   [3]  0.633333      0.1229175       ruleset1
#>   [4]  0.666667      0.0214385       ruleset1
#>   [5]  0.700000      0.2763628       ruleset1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### A note on primer design

The primer design component of the `mutateR` pipeline relies on the
`primer3-py` module accessed via the Python backend, and requires the
correct installation and activation of the `r-mutater` env.

This functionality applies two different primer design strategies
depending on the expected size of the deletion evoked by a given gRNA
pair:

- **Strategy A** is a simple, single primer pair approach used for when
  the wild type amplicon would be \< 3 kb. The forward and reverse
  primers are positioned slightly distal to the deletion sites,
  producing a large wild type amplicon and a small mutant one.

- **Strategy B** is a two-pair approach used for very large deletions
  (where the wild type amplicon size using Strategy A primers would be
  unreasonably large). The first pair flanks the deletion site, but will
  only produce an amplicon when the deletion is present (proximalising
  the primer binding sites). The second primer pair (nested within the
  deleted region) is used to distinguish wild type and heterozygous
  genotypes from homozygous ones.

The following parameters are the ones passed to `primer3-py`:

``` python
base_global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 65.0,
        'PRIMER_MAX_POLY_X': 4,       # Disallows homopolymers >4 nts
        'PRIMER_GC_CLAMP': 1,         # Weighted preference for 3' GC clamp
        'PRIMER_MAX_SELF_ANY': 8.00,
        'PRIMER_MAX_SELF_END': 3.00,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_NUM_RETURN': 1
    }
```

These are not user-facing (changing them requires modification of the
source code). However, the user can adjust the parameters
`primer_max_wt` (the wild type amplicon size cutoff \[in nucleotides\]
before Strategy B is evoked, default `3000`) and `primer_tm` (the target
primer T<sub>m</sub>, default `60.0`) in `run_mutateR()`.

In order to significantly speed up the primer design process, I have
applied vectorisation to batch-pass the designed primers between R and
the Python backend.

### Special cases

#### When the same exon is targeted twice

In certain cases, both gRNAs in a pair will target the same exon (e.g.,
if two phase-compatible exons are separated by a single exon, such that
the intervening exon is to be removed). In these cases, specific rules
need to be applied in order to prevent the generation of small deletions
that would likely cause the retention of a frameshifted exon remnant in
the mature mRNA, and a consequent PTC.

To account for these, I have applied the following filters during gRNA
pair assembly for pairs where both gRNAs target the same exon:

- Do not recommend gRNA pairs that would be expected to cause a deletion
  of \<50 bp (small deletions can be annoying to genotype, can cause
  frameshifted exon retention + PTC, and if gRNAs target loci too close
  to each other they can cause steric hindrance between Cas RNPs,
  thereby preventing the required simultaneous DSBs).

- Only recommend gRNA pairs that would result in removal of \>70% of the
  original exon mass **OR** if the residual exon fragment is \<50 nt.

  - We calculate the residual exon mass (REM) by subtracting the
    expected genomic deletion size from the exon length, then determine
    a deletion ratio by dividing the deletion size by the exon length.

  - Removal of a large portion of the exon will increase the likelihood
    of removing important exonic splicing enhancers (ESEs), which will
    lead to the exclusion of the remaining exon fragment.

  - Residual exon fragments \<50 bp in length are frequently excluded by
    the spliceosome during mature mRNA assembly (see [Hwang and Cohen,
    MCB
    (1997)](https://www.tandfonline.com/doi/abs/10.1128/MCB.17.12.7099)).

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
#> PAM distribution:
#> 
#> AGG CGG GGG TGG 
#>  72  45  68 102
#> Scoring gRNAs using model: ruleset1
#> Computing on‑target scores using ruleset1 model...
#> Scored 287 guides using ruleset1.
#> Assembling valid gRNA pairs for opn4.1 ...
#> Assembling gRNA pairs for exon‑flanking deletions...
#> Detected probability-based scores (e.g. ruleset1). Using default cutoff: 0.5
#> Single-exon/two-exon gene detected: constructing intragenic deletion pairs.
#> Detected intragenic assembly mode (≤2 exons).
#> Designing genotyping primers for 465 recommended pairs...
#> Preparing primer design batch requests...
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> Running Batch Primer3 (465 designs)...
#> Mapping results to dataframe...
#> Designed primers for 465/465 pairs.
#> Plotting exon phase compatibility and gRNA pairs...
#> Error in combn(exon_ranks, 2) : n < m
#> mutateR pipeline completed for opn4.1, finding 465 gRNA pairs.
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

cut_site_5p
</th>

<th style="text-align:right;">

cut_site_3p
</th>

<th style="text-align:left;">

seqnames_5p
</th>

<th style="text-align:right;">

genomic_deletion_size
</th>

<th style="text-align:left;">

transcript_id
</th>

<th style="text-align:left;">

recommended
</th>

<th style="text-align:left;">

priming_strategy
</th>

<th style="text-align:left;">

primer_ext_fwd
</th>

<th style="text-align:left;">

primer_ext_rev
</th>

<th style="text-align:left;">

primer_int_fwd
</th>

<th style="text-align:left;">

primer_int_rev
</th>

<th style="text-align:left;">

exp_wt_size
</th>

<th style="text-align:right;">

exp_mut_size
</th>

<th style="text-align:right;">

exp_int_size
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

GATCATGCCCACTACATCAT
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.6101578
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

32843416
</td>

<td style="text-align:right;">

32843447
</td>

<td style="text-align:left;">

chr2
</td>

<td style="text-align:right;">

31
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

CACCATGATCCTCCACACCC
</td>

<td style="text-align:left;">

AAACATGTTTCCAGCGGTGC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

168
</td>

<td style="text-align:right;">

137
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

GTTTATAGTGAACCTAGCTG
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.6622205
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

32843416
</td>

<td style="text-align:right;">

32843569
</td>

<td style="text-align:left;">

chr2
</td>

<td style="text-align:right;">

153
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

CACCATGATCCTCCACACCC
</td>

<td style="text-align:left;">

TCTCGAATTTCCTTCCCCGC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

656
</td>

<td style="text-align:right;">

503
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

TTCAAACAATCCGAGCAGCG
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.6464406
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

32843416
</td>

<td style="text-align:right;">

32844024
</td>

<td style="text-align:left;">

chr2
</td>

<td style="text-align:right;">

608
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

CACCATGATCCTCCACACCC
</td>

<td style="text-align:left;">

CTGACTGGTGAGGGTTGGAC
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1026
</td>

<td style="text-align:right;">

418
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

AGTGTCAGTTCCCTGACCTT
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.5869812
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

32843416
</td>

<td style="text-align:right;">

32844551
</td>

<td style="text-align:left;">

chr2
</td>

<td style="text-align:right;">

1135
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

CACCATGATCCTCCACACCC
</td>

<td style="text-align:left;">

GCGTTGCAGGCATGTATCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1654
</td>

<td style="text-align:right;">

519
</td>

<td style="text-align:right;">

NA
</td>

</tr>

<tr>

<td style="text-align:left;">

TCCACACCCTTTCCCCACCG
</td>

<td style="text-align:left;">

CTAGTGGGCAGAAATCTGAG
</td>

<td style="text-align:right;">

0.5227872
</td>

<td style="text-align:right;">

0.7434274
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

32843416
</td>

<td style="text-align:right;">

32844642
</td>

<td style="text-align:left;">

chr2
</td>

<td style="text-align:right;">

1226
</td>

<td style="text-align:left;">

ENSDART00000018501
</td>

<td style="text-align:left;">

TRUE
</td>

<td style="text-align:left;">

Flanking
</td>

<td style="text-align:left;">

CACCATGATCCTCCACACCC
</td>

<td style="text-align:left;">

GCGTTGCAGGCATGTATCAG
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

NA
</td>

<td style="text-align:left;">

1654
</td>

<td style="text-align:right;">

428
</td>

<td style="text-align:right;">

NA
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
mutant in [Douek et al., IJMS
(2021)](https://www.mdpi.com/1422-0067/22/11/5948 "An Engineered sgsh Mutant Zebrafish Recapitulates Molecular and Behavioural Pathobiology of Sanfilippo Syndrome A/MPS IIIA")).
`mutateR` factors in predicted terminal exon PTCs when scanning for
tolerated deletions, and flags these cases in the `run_mutateR` output
in the `$pairs` dataframe, under the `terminal_exon_case` column.

## Using the `mutateR_viewer`

`mutateR_viewer` is a Shiny app you can use to further streamline the
gRNA pair identification workflow. First, execute `run_mutateR` with
`interactive = TRUE` to produce a plotly object in the `$plot` level.
Pass this plot object directly to `mutateR_viewer` to access the Shiny
app.

![](images/clipboard-1401140530.gif)

Currently, you can use the `mutateR_viewer` Shiny app to export all (or
a selection of) gRNA pairs calculated for a selected exon pair either as
plain text (by copying to the clipboard, as in the above GIF) or as a
CSV file.

The app’s data table reports the exons in the pair, the two protospacer
sequences, their on-target scores, the (approximate) genomic- and
transcript-specific deletion sizes, and genotyping primer sequences.
Additional metadata from the `$pairs` data.frame not displayed in the
app is included in the exported CSV.

## Manual function execution

You may wish to call lower-level functions directly (e.g., for use in
other programmatic pipelines). Note, however, that many of these
functions are built to take outputs from other `mutateR` functions.
Refer to the manual pages for specific input requirements for each
function.

## To be implemented

- More nucleases (maybe from CasPEDIA or similar?)

- All `crisprScore` scoring methods – Work in progress

- Integration of off-target prediction analysis

- Cross-species sequence/domain conservation scoring/visualisation

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
#>  [3] tibble_3.3.0                        ggplot2_4.0.1                      
#>  [5] purrr_1.2.0                         dplyr_1.1.4                        
#>  [7] BSgenome.Hsapiens.UCSC.hg38_1.4.5   BSgenome_1.77.2                    
#>  [9] rtracklayer_1.69.1                  BiocIO_1.19.0                      
#> [11] Biostrings_2.77.2                   XVector_0.49.1                     
#> [13] GenomicRanges_1.61.5                GenomeInfoDb_1.45.12               
#> [15] Seqinfo_0.99.2                      IRanges_2.43.2                     
#> [17] S4Vectors_0.47.2                    BiocGenerics_0.55.1                
#> [19] generics_0.1.4                      mutateR_0.0.1                      
#> [21] reticulate_1.44.1                  
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.2.3                   bitops_1.0-9               
#>  [3] httr2_1.2.1                 biomaRt_2.65.16            
#>  [5] rlang_1.1.6                 magrittr_2.0.4             
#>  [7] matrixStats_1.5.0           compiler_4.5.1             
#>  [9] RSQLite_2.4.5               dir.expiry_1.17.0          
#> [11] png_0.1-8                   systemfonts_1.3.1          
#> [13] vctrs_0.6.5                 stringr_1.6.0              
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
#> [35] stringi_1.8.7               Rcpp_1.1.0                 
#> [37] SummarizedExperiment_1.39.2 knitr_1.50                 
#> [39] Matrix_1.7-4                tidyselect_1.2.1           
#> [41] rstudioapi_0.17.1           dichromat_2.0-0.1          
#> [43] abind_1.4-8                 yaml_2.3.11                
#> [45] codetools_0.2-20            curl_7.0.0                 
#> [47] lattice_0.22-7              Biobase_2.69.1             
#> [49] withr_3.0.2                 KEGGREST_1.49.1            
#> [51] S7_0.2.1                    evaluate_1.0.5             
#> [53] crisprScoreData_1.13.0      BiocFileCache_2.99.6       
#> [55] xml2_1.5.0                  ExperimentHub_2.99.5       
#> [57] pillar_1.11.1               BiocManager_1.30.27        
#> [59] filelock_1.0.3              MatrixGenerics_1.21.0      
#> [61] crisprScore_1.13.1          RCurl_1.98-1.17            
#> [63] BiocVersion_3.22.0          hms_1.1.4                  
#> [65] scales_1.4.0                glue_1.8.0                 
#> [67] tools_4.5.1                 AnnotationHub_3.99.6       
#> [69] GenomicAlignments_1.45.4    XML_3.99-0.20              
#> [71] grid_4.5.1                  AnnotationDbi_1.71.1       
#> [73] basilisk_1.21.5             restfulr_0.0.16            
#> [75] cli_3.6.5                   rappdirs_0.3.3             
#> [77] textshaping_1.0.4           kableExtra_1.4.0           
#> [79] viridisLite_0.4.2           S4Arrays_1.9.1             
#> [81] svglite_2.2.2               gtable_0.3.6               
#> [83] digest_0.6.39               SparseArray_1.9.1          
#> [85] rjson_0.2.23                farver_2.1.2               
#> [87] memoise_2.0.1               htmltools_0.5.9            
#> [89] lifecycle_1.0.4             httr_1.4.7                 
#> [91] bit64_4.6.0-1
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

Scott, D. A. & Zhang, F. Implications of human genetic variation in
CRISPR-based therapeutic genome editing. *Nat Med* **23**, 1095–1101
(2017). <https://doi.org/10.1038/nm.4377>

Yang, W., Zhu, J. K. & Jin, W. A catalog of gene editing sites and
genetic variations in editing sites in model organisms. *BMC Genomics*
**25**, 1153 (2024). <https://doi.org/10.1186/s12864-024-11073-9>

Kim, H. K. *et al.* Deep learning improves prediction of CRISPR-Cpf1
guide RNA activity. *Nat Biotechnol* **36**, 239–241 (2018).
<https://doi.org/10.1038/nbt.4061>

Kim, H. K. *et al.* SpCas9 activity prediction by DeepSpCas9, a deep
learning-based model with high generalization performance. *Sci Adv*
**5**, eaax9249 (2019). <https://doi.org/10.1126/sciadv.aax9249>

Hwang, D. Y. & Cohen, J. B. U1 small nuclear RNA-promoted exon selection
requires a minimal distance between the position of U1 binding and the
3’ splice site across the exon. *Mol Cell Biol* **17**, 7099–7107
(1997). <https://doi.org/10.1128/MCB.17.12.7099>

Douek, A. M. *et al.* An Engineered sgsh Mutant Zebrafish Recapitulates
Molecular and Behavioural Pathobiology of Sanfilippo Syndrome A/MPS
IIIA. *Int J Mol Sci* **22**, 5948 (2021).
<https://doi.org/10.3390/ijms22115948>
