# hilldiv3 manuscript data

This repository stages the Microflora Danica data panel for a worked example in
the **hilldiv3** manuscript. It contains a genome by sample abundance table,
sample and genome metadata, and a phylogenetic tree on the same 5,518 genomes.
Functional annotations and hilldiv3 analyses will be added in later stages.

## Bookdown publication

This repository follows the [alberdilabr](https://github.com/alberdilab/alberdilabr)
bookdown layout. `index.Rmd` is the landing page and defines the site output;
the numbered `.Rmd` files in the repository root are chapters, ordered by
`_bookdown.yml`. The current chapter documents the staged data. As analyses
are completed, add their chapters to `_bookdown.yml`. Shared setup, styling,
references, and generated figures live in `alberdilabr/`.

R package versions are pinned in `renv.lock`. After cloning, run
`renv::restore()` to install them. To build the site locally, run:

```sh
Rscript -e 'bookdown::render_book("index.Rmd")'
```

The site is generated in ignored `_site/`; HTML is not committed. Pushing to
`main` runs `.github/workflows/publish.yml`, which builds the site and deploys
it to GitHub Pages. Pages is configured to use GitHub Actions, and the
publication is available at <https://alberdilab.github.io/hilldiv3-ms/> after
the first successful deployment. When adding analysis dependencies, run
`renv::snapshot()` and commit the updated lockfile.

## Staged files

| File | Contents |
| --- | --- |
| `data/mag_relative_abundance.tsv.xz` | 5,518 representative MAGs (rows) × 360 samples (columns); first column is `genome_id`. Values are proportions of the published MFD MAG profile. |
| `data/sample_metadata.tsv` | One row per sample, in abundance-column order, including the eight habitat classes, MFD ontology, location, date, and selection audit. |
| `data/representative_metadata.tsv` | One row per 95% ANI species representative, in abundance-row order, including cluster, taxonomy, and genome-quality fields. |
| `data/representative_tree.nwk` | Newick tree pruned to exactly those 5,518 `genome_id` tips, with source branch lengths retained. |

The 360 samples are the geographically dispersed panel locked for the
[gifter Microflora Danica case study](https://github.com/alberdilab/gifter/tree/f04f7de0691c29874f7fbe17b915990b832e2a4d/manuscript/analysis/r10-mfd):
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

The representative and sample tables are byte-for-byte copies of the locked
panels in [gifter commit
`f04f7de`](https://github.com/alberdilab/gifter/tree/f04f7de0691c29874f7fbe17b915990b832e2a4d/manuscript/analysis/r10-mfd).
The source tree is
[`gtdb-shallow.tree`](https://github.com/cmc-aau/mfd_sr_mags/blob/957591d4104462ec6e52ab58d85a69798057e582/analysis/datasets/gtdb-shallow.tree)
from the MFD companion repository at commit
`957591d4104462ec6e52ab58d85a69798057e582` (SHA-256
`693976bea56093394704861c8229a26305a6aec21d8f60fde9debedc89acc112`).
The source tree includes other reference genomes; this repository removes
their tips and keeps the MFD representatives.

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
