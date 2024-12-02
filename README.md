## The Seroepidomiology of Non-Typhoidal Salmonella (NTS)

This project aims to explain the genomic findings of exploring NTS Salmonella epidemiology in Africa and the combination with ELISA results and the visualizations used in that project. This project involves several steps:

**1. Running bactopia**
- First create a fofn of samples via:
  ```
  bactopia prepare -p . > samples.txt
  ```
  (assuming the current directory has compressed paired-end FASTQs for the samples)

- Run bactopia via:
  ```
  bactopia --samples samples.txt -profile docker  --outdir bactopia_analysis
  ```

**2. Downstream analysis pipeline for Salmonella seroepidemiology**
- Library imports
  - Loads the required R libraries
- Data loading
  - Loads the required data e.g SISTR results in excel format
- Data preprocessing
  - Preprocess the data for downstream analysis e.g splitting into subspecies I & II
- Statistical analyses
  - Performs statistical analyses on the preprocessed data
- Visualization generation
  - Generates visualizations for the preprocessed data e.g histograms of serovars distributed in subspecies I & II