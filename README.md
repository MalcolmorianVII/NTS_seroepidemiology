## The Seroepidomiology of Non-Typhoidal Salmonella(NTS)
This project aims to explain the genomic findings of exploring NTS Salmonella epidemiology in Africa and the combination with ELISA results and the visualizations used in that project.This project involves several steps:
`1. ` Running bactopia 
    - First create a fofn of samples via
      - ```bactopia prepare -p . > samples.txt``` assuming the current directory has compressed paired end fastqs for the samples
    - Run bactopia via
      - ```bactopia --samples samples.txt -profile docker  --outdir bactopia_analysis```
`2.` Downstream analysis
    - Combine assembly-scan with SISTR results
    - Visualize the genomic findings in R
      - Distribution of serovars among subspecies I & II
      - Distribution of serogroups among subspecies I & II
    - Visualize ELISA results with genomic findings 
      - Visualize ELISA results with serovar classification
      - Visualize ELISA results with serogroup classification
      - Visualize ELISA results with unique antigenic profile classification