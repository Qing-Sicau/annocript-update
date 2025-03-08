# annocript-update
This is an updated python version of annocript

the original version was wrote by **Francesco Musacchia**

please cite: Annocript: a flexible pipeline for the annotation of transcriptomes able to identify putative long noncoding RNAs, doi: Bioinformatics, 2015 31(13):2199-201;

doi: 10.1093/bioinformatics/btv106

this is the AI automated generated project with only minor mamually curation, do not use under industry conditions.


General point:
- no downloading process, download all database files by your self.
- use diamond to substitute blastp
- no mysql

All scripts were generated using Grok.

# how to use:
``conda create -n annocript python=3.9``

``conda activate annocript``

``conda install -c bioconda diamond blast transdecoder pyarrow``

``install pandas dask duckdb pyfaidx biopython pyyaml``

git all content from this repo
``mkdir ./database``
### data preparation
- put all your downloaded databses files into this directory.
\#\#- from CPAT https://cpat.readthedocs.io/en/latest/, download Arabidopsis_hexamer.tsv and Arabidopsis_logit.RData, move into ./databases/
## annotate your fasta
- before formal run,
  ```bash
  python annocript.py
  ```
- then modify the parameters using vi or other tools as your favorate
- run the tool, using commands like:
```
python annocript.py --fasta plant_transcripts.fasta --threads 16 --do_blastx --do_blastn --do_rpstblastn --do_lnc_prediction --do_dna2pep --do_build_output --extract_stats
```
find your results in the out directory.
