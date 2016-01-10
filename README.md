# Tissue/organ specific gene expression data

`human-bodymap2.0-fpkms.bash` - a script to process Human BodyMap data, by Eric Minkel. [Source](https://gist.github.com/ericminikel/7533289). His post, ["Tissue-specific gene expression data based on Human BodyMap 2.0"](http://www.cureffi.org/2013/07/11/tissue-specific-gene-expression-data-based-on-human-bodymap-2-0/) provides details.

## `Gilad` folder

Supplementary data from Gilad Y, Mizrahi-Man O: [A reanalysis of mouse ENCODE comparative gene expression data](http://f1000research.com/articles/4-121/v1). F1000Research 2015, 121:1–32. [Download](https://zenodo.org/record/17606)

- `dataAndConfigurationFiles` - human/mouse gene annotation GTF files, GC content of non-overlapping exons, human/mouse one-to-one gene orthologs.
- `python_scripts` - scripts to generate some of the human/mouse gene annotation files
- `R_input_files` - Tissue-specific gene expression data as FPKMs (`Stanford_datasets_fpkmMat.txt`) and raw counts (`Stanford_datasets_rawCountsMat.txt`), 14,744 genes vs. 26 samples, no header/row names. `ortholog_GC_table.txt` - header plus human/mouse gene names and GC content, 4 columns. `Stanford_datasets.txt` - annotations of the 26 samples,  tissue, organism, and batch effect, 4 columns. Not corrected for batch.
- `R_scripts` - R scripts used to generate figures for the paper. `Supplementary_text_1/3.txt` - scripts to recreate published workflow. `Supplementary_text_2.txt` - accounting for batch effect.

`R_scripts/work_script_Gilad.R` - a short version of data processing pipeline. The depth-normalized, GC- and batch adjusted log2-transformed counts data is in `R_input_files/Stanford_dataset_normalized_counts_full.txt`


