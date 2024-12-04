### The Seroepidemiology of Non-Typhoidal Salmonella (NTS)

This project aims to integrate genomic analyses with serological data to investigate the epidemiology of non-typhoidal Salmonella (NTS) in Africa. By combining findings from whole-genome sequencing (WGS) with ELISA data and visualizations, the project provides insights into the seroepidemiology of NTS. The workflow involves the following steps:

---

#### **1. Running Bactopia**
Bactopia is used for processing raw WGS data. The steps include:

- **Prepare a file of filenames (fofn):**
  Generate a list of samples, assuming the current directory contains compressed paired-end FASTQ files:
  ```bash
  bactopia prepare -p . > samples.txt
  ```

- **Run the Bactopia pipeline:**
  Execute the analysis using the generated `samples.txt`:
  ```bash
  bactopia --samples samples.txt --profile docker --outdir bactopia_analysis
  ```

---

#### **2. Downstream Analysis Pipeline for Salmonella Seroepidemiology**
This pipeline integrates genomic and serological data to extract meaningful insights.

- **Library Imports:**
  Load the necessary R libraries for data manipulation and visualization.

- **Data Loading:**
  Import relevant datasets such as SISTR results (in Excel or other formats).

- **Data Preprocessing:**
  Prepare data for analysis, including tasks such as:
  - Splitting data into **subspecies I** and **subspecies II**.
  - Filtering and cleaning data.

- **Statistical Analyses:**
  Perform statistical analyses to explore trends and associations in the data.

- **Visualization Generation:**
  Create visual representations of the data, such as:
  - Histograms of serovar distributions across subspecies.
  - Other graphical summaries of genomic and serological findings.

---

This workflow ensures a comprehensive approach to understanding the interplay between genomic and serological data, providing valuable insights into NTS epidemiology in Africa.