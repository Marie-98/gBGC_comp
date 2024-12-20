# gBGC_comp

This Github project contains the scripts and files to perform all the pipelines and analyses for the project "Molecular re-adaptation : compensatory evolution following deleterious episodes of GC-biased gene conversion in rodents". 

Two main scripts contain the commands to be executed :
- pipeline_part1.sh : contains all the commands from reads to cleaned alignments
- pipeline_part2.sh : contains all the commands from cleaned alignments to final tables for analyses
In these scripts you have all the steps with all the commands used. Comments indicate which tools need to be installed for each command. At some points of the pipeline, there are comments recommending the suppression of some intermediates useless files that take a lot of memory. Some commands call for scripts : all these called scripts are available in the "scripts" folder. 
You will also need some in files which are all in the "in_files" folder, except for the reference exon sequences file, which is available here : https://zenodo.org/records/14534843?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImUxYmRiNTJhLTA3M2UtNGY2Yi05ZTc3LTlkZWU0MmM1OTIyOCIsImRhdGEiOnt9LCJyYW5kb20iOiIwOTM5NTNkNDBiZThkOTk4ZDI5ZTA1NDAwNWUzNjgzMyJ9.m9yPdEfreGPY9NHxMtik1vRGLHNNNLwXjsveXk08d7fqZL2IpX391lTXARbyNllB8Kj3664U4wmxFOW7EZUD-w. Comments in the pipeline explain where to put all these in files.

As our pipeline is quite long to run (~3 weeks), all the final files we used for our analyses are in the "final_files_for_analyses" folder, available here : https://zenodo.org/records/14534843?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImUxYmRiNTJhLTA3M2UtNGY2Yi05ZTc3LTlkZWU0MmM1OTIyOCIsImRhdGEiOnt9LCJyYW5kb20iOiIwOTM5NTNkNDBiZThkOTk4ZDI5ZTA1NDAwNWUzNjgzMyJ9.m9yPdEfreGPY9NHxMtik1vRGLHNNNLwXjsveXk08d7fqZL2IpX391lTXARbyNllB8Kj3664U4wmxFOW7EZUD-w.


