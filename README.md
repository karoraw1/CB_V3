# CB_V3

## OTU Processing Steps

### Dependencies:

1. Miniconda3 (for installing python and packages)
2. Python 3.7.3
3. numpy 1.16.3
4. pandas 0.24
5. pyyaml 5.1
6. biopython 1.73
7. R 3.5.1
8. dada2 1.11.1
9. FASTQC

### Process:

1. In `otu_scripts` there is a file called `config.yml`. Fill this out first.
2. Make sure paths are correct in the `*skeleton.sh` files in `utility_scripts`
3. Make sure your trimming parameters are set in the `TrimmingParams.tsv` file
4. Run `python prepSSnMakeJobs.py` to create batch files for each library
5. Run all the batch files that were created to trim your libraries