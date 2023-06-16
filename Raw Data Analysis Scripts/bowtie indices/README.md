# Bowtie Indices

Bowtie indices are the files that [bowtie](https://bowtie-bio.sourceforge.net/index.shtml) uses to align your fastq reads to the reference genome of your choosing. A single index consists of 6 files with the endings .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, rev.1.ebwt, rev.2.ebwt. They are [easy to make](https://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer)!  And all the indices that were used in the project can be found in this folder.  Here is a description of what each one is:

---
<br>
<details>
<summary>bsub_ecoli_capsules</summary>
<br>

This is a joint index for B subtilis (168, NCBI record: NC_000964.3), E coli (MG1655, NCBI record: NC_000913.2) and the capsule plasmid (pMP025).

By making a joint index I am able to simultaneously align my reads from the "mixing experiment" - where I mixed the capsule-containing E coli samples with B. sub, to all three genome references and discard any reads that do not uniquely map to one position in the genome (meaning we can disregard those reads that may align to rRNA or other components that remain similar between the two species.)

</details>
<br>
<details>
<summary>bsubtilis</summary>
<br>


This index is made from NCBI record NC_000964.3 for BS 168.

</details>
<br>
<details>
<summary>syn_ecoli</summary>
<br>

This is a joint index made from both E coli (MG1655, NCBI record: NC_000913.2) and the capsule plasmid (pMP025).

It allows us to map read originating from the capsule-containing E. coli strain used in these experiments.

</details>