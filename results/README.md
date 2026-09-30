# Results

Final isoform models and read counts behind manuscript Figs. 6–7 and Tables 2–4. File names
reflect the workflow steps that produced them.

```
results/
├── KLK4/                    KLK4 models without the KLKP1 fusion part ("nofusion")
│   ├── dRNA/
│   │   ├── all/             27 KLK4 models of the 30 final dRNA models, all 8 samples (Fig. 6A/B)
│   │   │   ├── *.counts.renamed.nofusion.v2.tsv   read counts per sample
│   │   │   ├── *predictProdORF.nofusion.v2.bed    model structure (BED12, modified chr19) + ORF tag
│   │   │   ├── fasta_models/                      sequence per model (hap1 if phased)
│   │   │   └── fasta_models_hap2/                 haplotype-2 sequences of the 12 phased models
│   │   ├── all-top10/       10 most abundant models
│   │   ├── 22/, 22-top10/   22Rv1 samples only
│   │   ├── C4/, C4-top10/   C4-2B samples only
│   │   └── phasing/phasing_results.v2.tsv        per-haplotype counts (step 09)
│   ├── cDNA/                same layout, 23 KLK4 models of the 37 final cDNA models (Fig. 7A/B)
│   └── comp_dRNA_cDNA/      models found by both methods / by one method only (Table 3)
├── KLK4T2/                  KLK4T2-like fusion models: 3 dRNA, 14 cDNA (Figs. 6C/D, 7C/D), same layout
├── KLK3/                    KLK3 models after step 06 (FLAIR quantification; not curated)
└── KLK2/                    KLK2 models after step 06 (FLAIR quantification; not curated)
```

* BED coordinates refer to the **modified chr19** (`chr19-cleaned-mod2.fa`, see step 00).
  Positions > chr19:50,906,558 are 67 bp larger than in GRCh38.
* Count columns: `22_0W 22_0 22_10W 22_10 C4_0W C4_0 C4_10W C4_10` (`W` = whole-cell
  lysate, no suffix = pooled polysome fractions).
* dRNA counts are FLAIR counts of the curated models. cDNA counts in `all/` are the
  alignment-based counts (`bamQuant`).
* KLK2/KLK3: the manuscript uses the KLK3 read counts (Tables 3 and 4). These models did not
  go through the curation of steps 07–09.
