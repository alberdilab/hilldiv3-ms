# hilldiv3 manuscript data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24001411.svg)](https://doi.org/10.5281/zenodo.24001411)

> **hilldiv3 is an R package for neutral, phylogenetic and functional
> Hill-number diversity analysis.** Source code, documentation and issue
> tracker: **[github.com/alberdilab/hilldiv3](https://github.com/alberdilab/hilldiv3)**
> · Package website: [alberdilab.github.io/hilldiv3](https://alberdilab.github.io/hilldiv3/)
> · This webbook: [alberdilab.github.io/hilldiv3-ms](https://alberdilab.github.io/hilldiv3-ms/)

This repository stages the Microflora Danica data panel for a worked example in
the **hilldiv3** manuscript. It contains a genome by sample abundance table,
sample and genome metadata, and a phylogenetic tree on the same 5,518 genomes.
The webbook includes a data overview, neutral and phylogenetic diversity
analyses, diversity profiles across orders, and pairwise beta diversity
comparisons. It also allocates Hill power-sum shares to taxonomic sets with
`hillshare()`. Its `hilldiv3` calls are highlighted separately from the plotting
and tabulation code. A focused `hillcontrib()` comparison attributes
agricultural field versus natural grassland soil dissimilarity to the same
taxonomic sets. The chapter also examines opposing Gaiellales and
Streptosporangiales patterns within a nearly stable Actinomycetota share.
A final chapter adds functional diversity from 17 function-level and 161
element-level traits derived from KEGG ortholog annotations, and tests whether
that order-level turnover changes what the community can do.

## Bookdown publication

The webbook is built from numbered R Markdown chapters in the repository root,
styled with the AlberdiLab GitBook theme in `style.css` and `_output.yml`, and
rendered into `docs/`. The chapters cover data preparation, a geographic and
genomic overview, initial analyses with `hilldiv3`, pairwise beta diversity,
taxonomic power-sum shares and dissimilarity contributions, and functional
diversity.

Install R packages `bookdown`, `rmarkdown`, `knitr`, `data.table`, `ape`,
`ggplot2`, `maps`, and `hilldiv3` 3.2.0 or later. The webbook in `docs/` was
rendered with `hilldiv3` 3.2.0 (source commit `3f2f98f`). Install it from
CRAN with `install.packages("hilldiv3")`, or from source with
`devtools::install_github("alberdilab/hilldiv3")`.
To render the webbook from the repository root, run:

```sh
Rscript scripts/render-webbook.R
```

The render script embeds the existing Fuse.js search dependency so the
`self_contained` pages do not depend on a CDN at viewing time.

Commit the updated `docs/` output with the R Markdown sources. GitHub Pages
serves `docs/` from `main` at
<https://alberdilab.github.io/hilldiv3-ms/>.

## Staged files

| File | Contents |
| --- | --- |
| `data/mag_relative_abundance.tsv.xz` | 5,518 representative MAGs (rows) × 360 samples (columns); first column is `genome_id`. Values are proportions of the published MFD MAG profile. |
| `data/sample_metadata.tsv` | One row per sample, in abundance-column order, including the eight habitat classes, MFD ontology, location, date, and selection audit. |
| `data/representative_metadata.tsv` | One row per 95% ANI species representative, in abundance-row order, including cluster, taxonomy, and genome-quality fields. |
| `data/representative_tree.nwk` | Newick tree pruned to exactly those 5,518 `genome_id` tips, with source branch lengths retained. |

The 360 samples are a geographically dispersed panel of
45 samples in each of eight habitat classes (agricultural field, natural
grassland, natural forest, urban greenspace, bog or fen soils; natural
freshwater and saltwater sediments; urban wastewater).

**These are relative abundances, not genome or read counts.** The source
Sylph profiles report percentages assigned to the MFD MAG collection. The
preparation script keeps strain (`t__`) rows, maps MAGs to their published
`secondary_cluster`, sums each 95% ANI cluster, and divides percentages by
100. It assigns the resulting cluster profile to the deposited representative
MAG. It does not renormalize profiles. A representative's abundance therefore
stands for its cluster and does not imply that all member MAGs have identical
gene content. These values describe the MFD MAG profile, not the fraction of
all DNA in each sample.

## Sources and attribution

The primary study is [Singleton and colleagues, *The Microflora Danica atlas
of Danish environmental microbiomes*](https://doi.org/10.1038/s41586-025-09794-2).
The [Microflora Danica data deposit](https://doi.org/10.5281/zenodo.17162544)
provides `Nitrifiers.tar.gz` (MD5
`9b65d9984cb41a929613758919ed2fcf`), which contains the published MAG
metadata and Sylph abundance table used here. Please cite both the study and
the deposit when using these data, and consult the deposit for its reuse terms.

## Licence and citation

This repository is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see `LICENSE`), the
same licence as the Microflora Danica data deposit it derives from. Cite the
release through `CITATION.cff` (GitHub's "Cite this repository") together with
the Microflora Danica study and deposit above. The archived releases are on Zenodo under the concept DOI
[10.5281/zenodo.24001411](https://doi.org/10.5281/zenodo.24001411), which always resolves to the latest version.

## Rebuild and check

The committed data are ready to use. To rebuild them from the pinned public
sources, install R packages `data.table` and `ape`, then run:

```sh
bash scripts/fetch-source-data.sh
Rscript scripts/prepare-data.R
Rscript scripts/check-data.R
```

The source archive is about 228 MB and its extracted abundance table is
compressed locally; downloaded inputs stay in the ignored `sources/`
directory. `prepare-data.R` also accepts three positional paths to an
already available MAG metadata table, abundance TSV.xz, and source tree, in
that order. To verify the committed snapshot byte-for-byte:

```sh
shasum -a 256 -c data/SHA256SUMS
```

The independent checker verifies the 360 sample IDs, 5,518 genome IDs, 45
samples per habitat class, abundance and metadata order, finite nonnegative
values, near-unit profile sums, and exact tree-tip identity. In this snapshot,
sample totals range from 0.999988 to 1.000015.
