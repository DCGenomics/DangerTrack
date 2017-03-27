# DangerTrack 🚨

A system for evaluating difficult to assess regions that uses SV calls and mappability to generate a genome-wide score.

---

Contents:

* Generated DangerTrack scores:
  * `dangertrack.bed.gz`
  * `dangertrack.bedgraph.gz`
* Scripts:
  * `prepare.sh` - step 1 - prepare, clean, and bin the data
  * `summarize.R` - step 2 - combine and summarize binned data
* Inputs:
  * `sv.1kg.bed.gz` - 1000 Genomes Project breakpoints
  * `sv.giab.bed.gz` - Genome in a Bottle breakpoints

---

originally part of https://github.com/NCBI-Hackathons/Structural_Variants_CSHL
