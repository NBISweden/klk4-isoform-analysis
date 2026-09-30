# Isoform naming conventions

Isoform names change three times during the workflow. This table lists every form you will
find in the scripts, parameter tables and result files.

| stage | example | meaning |
|---|---|---|
| `flair collapse` (step 04) | `ENST00000431178.2_ENSG00000167749.12`, `KLK4T2_…`, `dRNA-22_10-103_ENSG00000167749.12` | an annotated transcript, the KLK4T2 fusion transcript, or a representative read name (the sample name inside it has no meaning here), followed by the gene_id |
| after step 06 (`scripts/bin/rank_and_rename_models.py`) | `im1-a-R2344-P45.1_ENSG00000167749.12` | **i**soform **m**odel ranked 1 by read support; set **a** = corrected (annotated splice sites), **b** = inconsistent (novel splice sites); **R** = FLAIR read count; **P** = share of the reads of the set (set b: models with ≥ 10 reads), truncated to one decimal |
| step 08/09 (intermediate) | `im1-R2252` | final rank and read count of the FLAIR quantification of the curated models (step 08) |
| final / manuscript | `dim1-R2252`, `cim3-KLK4T2-R1393`, `dim8-ENST00000599865.5-0-R262` | **d**RNA or **c**DNA **i**soform **m**odel, rank, optional matching reference transcript, and total read count over all 8 samples (dRNA: FLAIR quantification; cDNA: alignment-based counts, step 08) |

The suffix `-0` after an Ensembl id (e.g. `ENST00000599865.5-0`) is added by FLAIR to
models assigned to that reference transcript (see the FLAIR documentation).

The ORF annotation from `flair predictProductivity` appears in the BED files of `results/` as
a suffix of the name, e.g. `dim1-R2252_PRO_ENSG00000167749.12`:

| tag | meaning |
|---|---|
| PRO | productive (ORF without premature termination codon) |
| PTC | premature termination codon (NMD candidate) |
| NGO | no start codon found |
| NST | start codon but no stop codon |

Haplotype sequences have the suffix `-hap1` / `-hap2` (`results/*/fasta_models*/`).

The mapping from step-06 names to the final names is in
[../config/dRNA_KLK4_step08_final_models.tsv](../config/dRNA_KLK4_step08_final_models.tsv) and
[../config/cDNA_KLK4_step08_final_models.tsv](../config/cDNA_KLK4_step08_final_models.tsv)
(column `final_name`).

## Sample names

| sample (manuscript Table 1) | material |
|---|---|
| `22_0`, `22_10`, `C4_0`, `C4_10` | pooled polysome fractions 12–21 |
| `22_0W`, `22_10W`, `C4_0W`, `C4_10W` | whole-cell lysate |

`22` = 22Rv1, `C4` = C4-2B, `0`/`10` = nM DHT.
