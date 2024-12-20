# gBGC_comp

This Github project contains the scripts and files to make all the pipeline and analyses for the project "Molecular re-adaptation : compensatory evolution following deleterious episodes of GC-biased gene conversion in rodents". 

Two main scripts contain the commands to run :
- pipeline_part1.sh : contains all the commands from reads to cleaned alignments
- pipeline_part2.sh : contains all the commands from cleaned alignments to final tables for analyses
In this scripts you have all the steps with all the commands used. Comments give the tools needed to be installed for each command. At some points of the pipeline, there are comments to recommend suppression of some intermediates useless files that take a lot of memory. Some commands call for scripts. All these called scripts are available in the "scripts" folder. 
You also need some in files that are all in the "in_files" folder, except the file of the sequences of reference exons that is available here : https://zenodo.org/records/14534843?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImUxYmRiNTJhLTA3M2UtNGY2Yi05ZTc3LTlkZWU0MmM1OTIyOCIsImRhdGEiOnt9LCJyYW5kb20iOiIwOTM5NTNkNDBiZThkOTk4ZDI5ZTA1NDAwNWUzNjgzMyJ9.m9yPdEfreGPY9NHxMtik1vRGLHNNNLwXjsveXk08d7fqZL2IpX391lTXARbyNllB8Kj3664U4wmxFOW7EZUD-w. Comments in the pipeline explained where to put all these in files.

As our pipeline is very long to execute (~3 weeks), all the final files that we used to make our analyses are in the "final_files_for_analyses" folder available here : https://zenodo.org/records/14534843?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImUxYmRiNTJhLTA3M2UtNGY2Yi05ZTc3LTlkZWU0MmM1OTIyOCIsImRhdGEiOnt9LCJyYW5kb20iOiIwOTM5NTNkNDBiZThkOTk4ZDI5ZTA1NDAwNWUzNjgzMyJ9.m9yPdEfreGPY9NHxMtik1vRGLHNNNLwXjsveXk08d7fqZL2IpX391lTXARbyNllB8Kj3664U4wmxFOW7EZUD-w.


