**Project: Closing genome of *Tenacibaculum mesophilum* strain ECR**  
**Authors: Rebecca Mickol and Artemis Louyakis**  
**Date: 2020.01.07**  


Combined the two illumina runs into single fwd and rev files. Ran unicycler assembly optimizer (spades for illumina / miniasm+racon for long reads) in a hybrid assembly first with the unassembled minion data, then with the assembled data. Additionally, ran unicycler_polish using the minion assembly and polished with the illumina reads. The latter resulted in the best overall assembly, so far. CheckM was used to test the assembly with 100% completion using the family marker set and 92.67% completion with the genus marker set. The low genus score is likely because only two genomes were used to generate the markers and it's likley that the species are more divergent than those two genomes can account for. 

Program versions:
unicycler v0.4.7
Racon v0.5.0
SPAdes v3.13.0
samtools v1.9
minimap2 v2.17
bowtie2 v2.3.5.1
freebayes v1.1.0
PBSuite v15.8.24
pilon v1.22
blast v2.7.1

#### Housekeeping:
```bash
cat illumina_run1/*R1* illumina_run2/*R1* > ecto_r1.fastq
cat illumina_run1/*R2* illumina_run2/*R2* > ecto_r2.fastq

ls
# illumina_run2  ecto_r1.fastq  ecto_r2.fastq  illumina_run1  minion

ls minion/
# ecto.contigs.fasta  ecto.unassembled.fasta

### calculate N50
grep "^>" minion/ecto.unassembled.fasta | awk -F" " '{print $2}' | cut -d"=" -f2 | ./n50.py 
# 5623
```

#### Unicycler: hybrid assemblies with unassembleda and assembled long reads
```bash
unicycler -1 ecto_r1.fastq -2 ecto_r2.fastq -l minion/ecto.unassembled.fasta \
  -o scratch -t 10 --mode bold
# 47 contigs; longest 971,344 bp

unicycler -1 ecto_r1.fastq -2 ecto_r2.fastq -l minion/ecto.contigs.fasta \
  -o scratch2 -t 10 --mode bold
# 28 contigs; longest 1,199,730

unicycler -1 ecto_r1.fastq -2 ecto_r2.fastq -l minion/ecto.unassembled.fasta \
  --existing_long_read_assembly minion/ecto.contigs.fasta \
  -o scratch3 -t 10 --mode bold
# 47 contigs; longest 971,344 bp

unicycler -1 illumina_run2/LK_RLM_AGGCAGAA-GTAAGGAG_L001_R1_001.fastq -2 illumina_run2/LK_RLM_AGGCAGAA-GTAAGGAG_L001_R2_001.fastq -l minion/ecto.unassembled.fasta -o illumina_2_unass -t 10 --mode bold
# 43 contigs; longest 1,739,238 bp

unicycler -1 illumina_run1/Sample-16_S16_L001_R1_001.fastq -2 illumina_run1/Sample-16_S16_L001_R2_001.fastq -l minion/ecto.unassembled.fasta -o illumina_1_unass -t 10 --mode bold
# 104 contigs; longest 210,188 bp

unicycler -1 ../illumina_run1/Sample-16_S16_L001_R1_001.fastq -2 ../illumina_run1/Sample-16_S16_L001_R2_001.fastq -l minion/ecto.contigs.fasta -o illumina_1 -t 10 --mode bold
# 49 contigs; longest 1,457,567 bp

unicycler -1 illumina_run2/LK_RLM_AGGCAGAA-GTAAGGAG_L001_R1_001.fastq -2 illumina_run2/LK_RLM_AGGCAGAA-GTAAGGAG_L001_R2_001.fastq -l minion/ecto.contigs.fasta -o illumina_2 -t 10 --mode bold
# 24 contigs; longest 1,964,511 bp

unicycler -1 trim_ecto_r1.fq -2 trim_ecto_r2.fq -s trim_ecto_s.fq -l minion/ecto.contigs.fasta -o uni_trim/ --mode bold
# 24 contigs; longest 1,964,514 bp

## output available upon request
```

#### Unicycler: polished illumina assembly with illumina reads
```bash
### this assembly used contigs from a spades assembly of the illumina_run2 data polished with all the trimmed illumina reads
cd polish_spades
unicycler_polish -a ../new_raw/spades_16/scaffolds.fasta \
  -1 ../trim_ecto_r1.fq -2 ../trim_ecto_r2.fq --threads 10 \
  --ale ~/scripts/ALE/src/ALE \
  --samtools ~/scripts/miniconda/bin/samtools \
  --bowtie2 ~/scripts/miniconda/bin/bowtie2 \
  --pilon ~/scripts/miniconda/share/pilon-1.23-2/pilon-1.23.jar 
# 388 contigs; lengths 

### this assembly used spades assembly of all the illumina data polished with the trimmed illumina reads
cd polish_spades18
unicycler_polish -a ../2019_ecto_PE_SPAdes.contigs.FASTA/2019_ecto_PE_SPAdes.contigs.fa \
  -1 ../trim_ecto_r1.fq -2 ../trim_ecto_r2.fq --threads 10 \
  --ale ~/scripts/ALE/src/ALE \
  --samtools ~/scripts/miniconda/bin/samtools \
  --bowtie2 ~/scripts/miniconda/bin/bowtie2 \
  --pilon ~/scripts/miniconda/share/pilon-1.23-2/pilon-1.23.jar 
# 18 contigs (no change)

### Ran each with unassembled long reads to polish with and did not improve any of the assemblies; it's unclear if it even attempted the polishing - finished in a second
```

#### Unicycler: polished minion assembly with illumina reads
```bash
cd polish
unicycler_polish -a minion/ecto.contigs.fasta \
  -1 ecto_r1.fastq -2 ecto_r2.fastq --threads 10 \
  --ale ~/scripts/ALE/src/ALE \
  --samtools ~/scripts/miniconda/bin/samtools \
  --bowtie2 ~/scripts/miniconda/bin/bowtie2 \
  --pilon ~/scripts/miniconda/share/pilon-1.23-2/pilon-1.23.jar 
# 2 contigs; lengths 2083114 and 1322303
## described and uploaded version
## raw reads available on NCBI; trimmed illumina reads and assembled minion contigs are in data folder on github
```

```bash
ls ./bins/
# 014_final_polish.fasta
checkm taxon_list --rank family
checkm taxon_set family Flavobacteriaceae flavobacteriaceae.ms
checkm analyze flavobacteriaceae.ms ./bins ./out -x fasta -t 10
checkm qa -t 10 flavobacteriaceae.ms out/
```

*******************************************************************************
 [CheckM - qa] Tabulating genome statistics.
*******************************************************************************

  Calculating AAI between multi-copy marker genes.

  Reading HMM info from file.
  Parsing HMM hits to marker genes:
    Finished parsing hits for 1 of 1 (100.00%) bins.


| Bin Id                | Marker lineage     | #genomes  | #markers  | #marker sets  | 0   | 1   | 2  | 3  | 4  | 5+  | Completeness  | Contamination  | Strain heterogeneity  |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  014_final_polish  | Flavobacteriaceae (4)     | 115        | 457          | 309       | 0  | 455  | 2  | 0  | 0  | 0      | 100.00          | 0.32              | 0.00 |

```bash
checkm taxon_list --rank genus
checkm taxon_set genus Tenacibaculum tenacibaculum.ms
checkm analyze tenacibaculum.ms ./bins ./outg -x fasta -t 10
checkm qa -t 10 tenacibaculum.ms outg/
```
*******************************************************************************
 [CheckM - qa] Tabulating genome statistics.
*******************************************************************************

  Calculating AAI between multi-copy marker genes.

  Reading HMM info from file.
  Parsing HMM hits to marker genes:
    Finished parsing hits for 1 of 1 (100.00%) bins.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------
| Bin Id                | Marker lineage     | #genomes  | #markers  | #marker sets  | 0   | 1   | 2  | 3  | 4  | 5+  | Completeness  | Contamination  | Strain heterogeneity  |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|  014_final_polish   |Tenacibaculum (5)       |2          |1057          |272        |33   |985   |34   |5   |0   |0       |92.67            |4.27               |6.12 |


#### Align illumina reads to polished assembly to calculate alignment rate, then extract unaligned reads and convert to fastq format
```bash
bowtie2 -q --threads 10 -x polish/014_final_polish.fasta -1 ecto_r1.fastq -2 ecto_r2.fastq |\
  samtools view -b - |\
  samtools sort -@ 10 -m 5G -o coverage/014_final_polish.bam
samtools view -f4 coverage/014_final_polish.bam > coverage/unaligned_014_final_polish.sam

#cat unaligned_014_final_polish.sam | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > unaligned_014_final_polish_r1.fastq
#cat unaligned_014_final_polish.sam | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > unaligned_014_final_polish_r2.fastq
cat coverage/unaligned_014_final_polish.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > coverage/unaligned_014_final_polish.fq

```
```
bowtie2 alignment output:
601370 reads; of these:
  601370 (100.00%) were paired; of these:
    62088 (10.32%) aligned concordantly 0 times
    523791 (87.10%) aligned concordantly exactly 1 time
    15491 (2.58%) aligned concordantly >1 times
    ----
    62088 pairs aligned concordantly 0 times; of these:
      33374 (53.75%) aligned discordantly 1 time
    ----
    28714 pairs aligned 0 times concordantly or discordantly; of these:
      57428 mates make up the pairs; of these:
        35067 (61.06%) aligned 0 times
        20407 (35.53%) aligned exactly 1 time
        1954 (3.40%) aligned >1 times
97.08% overall alignment rate
```

```bash
bowtie2 -q --threads 10 -x polish/014_final_polish.fasta -1 trim_ecto_r1.fq -2 trim_ecto_r2.fq -U trim_ecto_s.fq |\
 samtools view -b - |\
 samtools sort -@ 10 -m 5G -o coverage/014_final_polish_trim.bam
samtools view -f4 coverage/014_final_polish_trim.bam > coverage/unaligned_014_final_polish_trim.sam
#cat unaligned_014_final_polish_trim.sam | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > unaligned_014_final_polish_trim_r1.fastq
#cat unaligned_014_final_polish_trim.sam | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > unaligned_014_final_polish_trim_r2.fastq
cat coverage/unaligned_014_final_polish_trim.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > coverage/unaligned_014_final_polish_trim.fq
```
```
bowtie2 alignment output:
601273 reads; of these:
  598163 (99.48%) were paired; of these:
    79723 (13.33%) aligned concordantly 0 times
    502831 (84.06%) aligned concordantly exactly 1 time
    15609 (2.61%) aligned concordantly >1 times
    ----
    79723 pairs aligned concordantly 0 times; of these:
      54010 (67.75%) aligned discordantly 1 time
    ----
    25713 pairs aligned 0 times concordantly or discordantly; of these:
      51426 mates make up the pairs; of these:
        30365 (59.05%) aligned 0 times
        18707 (36.38%) aligned exactly 1 time
        2354 (4.58%) aligned >1 times
  3110 (0.52%) were unpaired; of these:
    215 (6.91%) aligned 0 times
    2789 (89.68%) aligned exactly 1 time
    106 (3.41%) aligned >1 times
97.45% overall alignment rate
```

