# Silverside progeny relatedness test


## 

## Goal:

Áki’s candidate de novo SNP lists don’t make any sense when we’re
looking at the bam files in IGViewer, so I want to check if samples IDs
may have been scrambled and our indication of who is related how may be
mixed up

## Merge the two family VCFs into one 10-individual VCF

Áki has shared lists of separate vcfs of called sites with the bwa-GATK
pipeline for the two families.

`/fs/cbsubscb16/storage/silverside_progeny/ManualCheck/JP_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz`

`/fs/cbsubscb16/storage/silverside_progeny/ManualCheck/PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz`

I want to first create a vcf of shared SNPs (the intersection (SNPs
called in both), not union). I didn’t have write privileges to Áki’s
folder, so I’ve copied them to here:
`/workdir/nina_july2026/silverside_progeny`

``` bash

cd /workdir/nina_july2026/silverside_progeny/vcf

# First index bam files
tabix -p vcf JP_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz
tabix -p vcf PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz

bcftools merge \
  JP_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz \
  PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz \
  -Oz -o all10.merged.vcf.gz



tabix -p vcf all10.merged.vcf.gz
```

Take the intersection

``` bash

bcftools isec \
    -n=2 \
    -p shared_sites \
    JP_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz \
    PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz
```

This creates `shared_sites/sites.txt`

Keep only those sites shared between the families (the intersection) in
the merged VCF

``` bash

bcftools view \
    -R shared_sites/sites.txt \
    all10.merged.vcf.gz \
    -Oz \
    -o all10.shared.vcf.gz

tabix -p vcf all10.shared.vcf.gz
```

Check sample names

``` bash
bcftools query -l all10.shared.vcf.gz
```

## Keep only shared, clean, informative biallelic SNPs

This keeps SNPs that are:

- biallelic SNPs

- PASS variants

- no missing genotypes across all 10 individuals

- MAF \> 0.1 \[to remove noise from rare alleles - we’re looking at
  Mendelian relationships here, not detecting de novo mutations

I think Áki had already filtered for most of these criteria, so this is
just removing low-frequency variants

``` bash

bcftools view \
    -m2 -M2 \
    -v snps \
    -f PASS \
    -i 'F_MISSING==0 && MAF>0.1' \
    all10.shared.vcf.gz \
    -Oz \
    -o all10.shared.clean.snps.vcf.gz

tabix -p vcf all10.shared.clean.snps.vcf.gz
```

Check how many SNPs remain:

``` bash

echo "Before filtering JP only:"
bcftools view -H JP_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz | wc -l

echo "Before filtering PJ only:"
bcftools view -H PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz | wc -l

echo "Shared, before filtering:"
bcftools view -H all10.shared.vcf.gz | wc -l

echo "Shared, after filtering:"
bcftools view -H all10.shared.clean.snps.vcf.gz | wc -l
```

Before filtering JP only: 15217440

Before filtering PJ only: 15183988

Shared, before filtering: 14864656

Shared, after filtering: 8314322

## Convert to PLINK

I want to add SNP positions as the name for each site so I can track
where

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --vcf all10.shared.clean.snps.vcf.gz \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out all10_clean
```

## Edit the `.fam` file (created by plink in the conversion)

The format is

FamilyID IndividualID FatherID MotherID Sex Phenotype

I look at `nano all10.fam` I manually edited and changed this file to
contain

``` bash
cat all10_clean.fam 

JP LMJF3 0 0 2 -9
JP LMPM3 0 0 1 -9
JP LMJP001 LMPM3 LMJF3 0 -9
JP LMJP004 LMPM3 LMJF3 0 -9
JP LMJP005 LMPM3 LMJF3 0 -9
PJ LMPF2 0 0 2 -9
PJ LMJM2 0 0 1 -9
PJ LMPJ001 LMJM2 LMPF2 0 -9
PJ LMPJ002 LMJM2 LMPF2 0 -9
PJ LMPH003 LMJM2 LMPF2 0 -9
```

## Run Mendelian error analysis

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile all10_clean \
  --allow-extra-chr \
  --mendel \
  --out all10_clean_mendel
```

Useful outputs

``` bash

cat all10_clean_mendel.imendel

 FID     IID   N
  JP   LMPM3 1475016
  JP   LMJF3 2181533 
  JP LMJP001 380766
  JP LMJP004 386339
  JP LMJP005 1617236
  PJ   LMJM2 44286
  PJ   LMPF2 38861 
  PJ LMPJ001 23998
  PJ LMPJ002 23220
  PJ LMPH003 22445



cat all10_clean_mendel.fmendel
```

PLINK’s `.imendel` counts are **per individual involvement in Mendel
errors**, not “errors discovered independently in that person.”

A Mendel error is detected at the **trio level**: child + father +
mother. But when PLINK summarizes by individual, it increments the count
for **everyone in the trio involved in that inconsistent genotype**.

Interpretation of your output:

    JP: catastrophic inconsistency
    PJ: much lower, likely genotype noise

For PJ, the parents have higher counts than each child because each
parent is involved in **three offspring trios**. A bad parental genotype
can create Mendel errors with multiple offspring, so parents naturally
accumulate more errors.

For PJ, the pattern is compatible with a correct family plus some noisy
genotypes. It does **not** mean the parents themselves “failed Mendelian
segregation” in isolation.

The number of Mendelian errors for JP is so high that the provided
pedigree can’t be correct

## Estimate relatedness among all 10 individuals

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile all10_clean_ids \
  --allow-extra-chr \
  --genome \
  --out all10_clean_relatedness
```

Look at the output

``` bash
 FID1   IID1 FID2   IID2 RT    EZ      Z0      Z1      Z2  PI_HAT PHE       DST     PPC   RATIO
  JP  LMJF3  JP  LMPM3 OT     0  0.0000  1.0000  0.0000  0.5000  -1  0.783086  1.0000 115.6250
  JP  LMJF3  JP LMJP001 PO   0.5  0.0147  0.9853  0.0000  0.4927  -1  0.793101  1.0000 115.6250
  JP  LMJF3  JP LMJP004 PO   0.5  0.0000  1.0000  0.0000  0.5000  -1  0.784774  1.0000 154.3333
  JP  LMJF3  JP LMJP005 PO   0.5  1.0000  0.0000  0.0000  0.0000  -1  0.642896  0.0019  1.6459
  JP  LMJF3  PJ  LMPF2 UN    NA  0.9118  0.0584  0.0298  0.0590  -1  0.685600  0.9449  2.2396
  JP  LMJF3  PJ  LMJM2 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.642958  0.0000  1.5148
  JP  LMJF3  PJ LMPJ001 UN    NA  0.8391  0.1609  0.0000  0.0805  -1  0.637475  0.0478  1.7851
  JP  LMJF3  PJ LMPJ002 UN    NA  0.7864  0.2136  0.0000  0.1068  -1  0.647188  0.1237  1.8476
  JP  LMJF3  PJ LMPH003 UN    NA  0.8758  0.1242  0.0000  0.0621  -1  0.646253  0.7563  2.0997
  JP  LMPM3  JP LMJP001 PO   0.5  0.1402  0.4230  0.4368  0.6483  -1  0.867775  1.0000 19.3261
  JP  LMPM3  JP LMJP004 PO   0.5  0.1783  0.4495  0.3722  0.5969  -1  0.849763  1.0000 16.3148
  JP  LMPM3  JP LMJP005 PO   0.5  0.0141  0.9812  0.0047  0.4953  -1  0.796540  1.0000 115.6250
  JP  LMPM3  PJ  LMPF2 UN    NA  0.7508  0.2492  0.0000  0.1246  -1  0.637295  0.8227  2.1342
  JP  LMPM3  PJ  LMJM2 UN    NA  0.7715  0.2285  0.0000  0.1142  -1  0.636055  0.7415  2.0927
  JP  LMPM3  PJ LMPJ001 UN    NA  0.6952  0.3048  0.0000  0.1524  -1  0.624595  0.9449  2.2396
  JP  LMPM3  PJ LMPJ002 UN    NA  0.6924  0.3076  0.0000  0.1538  -1  0.632609  0.9929  2.3841
  JP  LMPM3  PJ LMPH003 UN    NA  0.7165  0.2835  0.0000  0.1417  -1  0.633441  0.9942  2.3964
  JP LMJP001  JP LMJP004 FS   0.5  0.2298  0.4370  0.3331  0.5517  -1  0.835199  1.0000 10.2651
  JP LMJP001  JP LMJP005 FS   0.5  0.0000  1.0000  0.0000  0.5000  -1  0.787264  1.0000 132.2857
  JP LMJP001  PJ  LMPF2 UN    NA  0.7590  0.2410  0.0000  0.1205  -1  0.640739  0.9902  2.3633
  JP LMJP001  PJ  LMJM2 UN    NA  0.7868  0.2132  0.0000  0.1066  -1  0.640215  0.9065  2.1952
  JP LMJP001  PJ LMPJ001 UN    NA  0.7059  0.2941  0.0000  0.1471  -1  0.628052  0.9102  2.1986
  JP LMJP001  PJ LMPJ002 UN    NA  0.7013  0.2987  0.0000  0.1493  -1  0.636216  0.8984  2.1877
  JP LMJP001  PJ LMPH003 UN    NA  0.7236  0.2764  0.0000  0.1382  -1  0.637779  0.9999  2.6484
  JP LMJP004  JP LMJP005 FS   0.5  0.0152  0.9848  0.0000  0.4924  -1  0.794955  1.0000 132.2857
  JP LMJP004  PJ  LMPF2 UN    NA  0.7800  0.2200  0.0000  0.1100  -1  0.646373  0.8040  2.1237
  JP LMJP004  PJ  LMJM2 UN    NA  0.8027  0.1973  0.0000  0.0986  -1  0.645871  0.4266  1.9745
  JP LMJP004  PJ LMPJ001 UN    NA  0.7206  0.2794  0.0000  0.1397  -1  0.633801  0.7340  2.0894
  JP LMJP004  PJ LMPJ002 UN    NA  0.7206  0.2794  0.0000  0.1397  -1  0.641760  0.8984  2.1877
  JP LMJP004  PJ LMPH003 UN    NA  0.7414  0.2586  0.0000  0.1293  -1  0.643384  0.9992  2.5113
  JP LMJP005  PJ  LMPF2 UN    NA  1.0000  0.0000  0.0000  0.0000  -1  0.646038  0.0001  1.5675
  JP LMJP005  PJ  LMJM2 UN    NA  0.9272  0.0177  0.0551  0.0639  -1  0.688729  0.4266  1.9745
  JP LMJP005  PJ LMPJ001 UN    NA  0.8033  0.1967  0.0000  0.0984  -1  0.643259  0.0186  1.7361
  JP LMJP005  PJ LMPJ002 UN    NA  0.8112  0.1888  0.0000  0.0944  -1  0.647562  0.5092  2.0032
  JP LMJP005  PJ LMPH003 UN    NA  0.8276  0.1724  0.0000  0.0862  -1  0.653434  0.0413  1.7768
  PJ  LMPF2  PJ  LMJM2 OT     0  1.0000  0.0000  0.0000  0.0000  -1  0.645413  0.0045  1.6762
  PJ  LMPF2  PJ LMPJ001 PO   0.5  0.0145  0.9855  0.0000  0.4928  -1  0.794068  1.0000 185.8000
  PJ  LMPF2  PJ LMPJ002 PO   0.5  0.0000  1.0000  0.0000  0.5000  -1  0.788272  1.0000 933.0000
  PJ  LMPF2  PJ LMPH003 PO   0.5  0.0000  1.0000  0.0000  0.5000  -1  0.787439  1.0000 185.6000
  PJ  LMJM2  PJ LMPJ001 PO   0.5  0.0000  1.0000  0.0000  0.5000  -1  0.786183  1.0000 102.6667
  PJ  LMJM2  PJ LMPJ002 PO   0.5  0.0168  0.9832  0.0000  0.4916  -1  0.792093  1.0000 92.5000
  PJ  LMJM2  PJ LMPH003 PO   0.5  0.0156  0.9844  0.0000  0.4922  -1  0.793013  1.0000 132.4286
  PJ LMPJ001  PJ LMPJ002 FS   0.5  0.2352  0.4917  0.2731  0.5189  -1  0.822340  1.0000 10.1310
  PJ LMPJ001  PJ LMPH003 FS   0.5  0.2588  0.5192  0.2220  0.4816  -1  0.808949  1.0000  9.7471
  PJ LMPJ002  PJ LMPH003 FS   0.5  0.3241  0.4852  0.1907  0.4333  -1  0.794172  1.0000  8.5408
```

The output reports pairwise estimates of identity-by-descent (IBD)
between every pair of individuals. The most informative columns are:

| Column | Interpretation |
|----|----|
| **IID1, IID2** | The two individuals being compared. |
| **PI_HAT** | Estimated proportion of the genome shared identical-by-descent. This is the primary statistic used to infer relatedness. |

Typical PI_HAT values are:

| Relationship                            | Expected PI_HAT            |
|-----------------------------------------|----------------------------|
| Same individual                         | ~1.0                       |
| Parent-offspring                        | ~0.5                       |
| Full siblings                           | ~0.5 (with some variation) |
| Half siblings / grandparent / avuncular | ~0.25                      |
| Unrelated                               | ~0                         |

## Interpretation of the PJ family

The PJ family shows the expected pattern for a correctly assigned
pedigree.

**Parent-offspring relationships:** All three offspring have PI_HAT
values of approximately 0.49–0.50 with both parents:

| Pair            | PI_HAT |
|-----------------|--------|
| LMPF2 – LMPJ001 | 0.493  |
| LMPF2 – LMPJ002 | 0.500  |
| LMPF2 – LMPH003 | 0.500  |
| LMJM2 – LMPJ001 | 0.500  |
| LMJM2 – LMPJ002 | 0.492  |
| LMJM2 – LMPH003 | 0.492  |

These values are exactly what is expected for true parent-offspring
relationships.

**Sibling relationships:** The three offspring are also related to one
another as expected:

Full siblings are expected to have PI_HAT ≈ 0.5, although individual
sibling pairs naturally vary because each sibling inherits a different
combination of parental chromosomes. These values therefore support the
expected full-sibling relationships.

**Between-family comparisons:** All comparisons between JP and PJ
individuals have PI_HAT values close to zero (approximately 0–0.15),
consistent with the two families being unrelated.

## Interpretation of the JP family

In contrast, the JP family is inconsistent with the recorded pedigree.

The expected pedigree is:

    LMJF3 + LMPM3
            ↓
    LMJP001, LMJP004, LMJP005

However, the observed relatedness does not support this.

### Evidence

LMJF3 has parent-offspring relatedness with LMJP001 and LMJP004:

| Pair            | PI_HAT |
|-----------------|--------|
| LMJF3 – LMJP001 | 0.493  |
| LMJF3 – LMJP004 | 0.500  |

but is essentially unrelated to LMJP005:

| Pair            | PI_HAT |
|-----------------|--------|
| LMJF3 – LMJP005 | 0.000  |

Similarly, LMPM3 has an expected parent-offspring relationship with
LMJP005:

| Pair            | PI_HAT |
|-----------------|--------|
| LMPM3 – LMJP005 | 0.495  |

but substantially higher-than-expected relatedness to LMJP001 and
LMJP004:

| Pair            | PI_HAT |
|-----------------|--------|
| LMPM3 – LMJP001 | 0.648  |
| LMPM3 – LMJP004 | 0.597  |

Finally, LMJF3 and LMPM3 themselves show a first-degree relationship:

| Pair          | PI_HAT |
|---------------|--------|
| LMJF3 – LMPM3 | 0.500  |

whereas the recorded pedigree indicates that they should be unrelated
parents.

### CONCLUSION

**It seems that the most likely explanation is that LMJP005 is the
father and LMPM3 is an offspring.**

With that swap, the key PI_HAT values make sense:

    LMJF3 – LMJP005   ≈ 0.00   unrelated parents
    LMJF3 – LMPM3     ≈ 0.50   parent-offspring
    LMJF3 – LMJP001   ≈ 0.49   parent-offspring
    LMJF3 – LMJP004   ≈ 0.50   parent-offspring

    LMJP005 – LMPM3   ≈ 0.50   parent-offspring
    LMJP005 – LMJP001 ≈ 0.50   parent-offspring
    LMJP005 – LMJP004 ≈ 0.49   parent-offspring

The only slight oddity is that `LMPM3` is unusually high with `LMJP001`
and `LMJP004` (~0.60–0.65), but that is much less problematic than the
original pedigree and could reflect related parents, marker
ascertainment, or genotype artifacts.

ChatGPT says it’s very unlikely to see such high pi_hat estimates with
millions of SNPs (~8 million) in offspring of unrelated, non-inbred
parents. So it’s suggesting that there might have been some mixing of
reads from different individuals when the bam files from the different
sequencing runs were merged. We need to investigate this.

## Mendelian analysis with updated pedigree

First make copy of the original .fam file and relatedness output

``` bash

cp all10_clean.fam all10_clean_original.fam
cp all10_clean_relatedness.genome all10_clean_original_relatedness.genome


# Then manually edit the .fam file to this (the sample order has to be the same)

JP LMJF3 0 0 2 -9
JP LMPM3 LMJP005 LMJF3 0 -9
JP LMJP001 LMJP005 LMJF3 0 -9
JP LMJP004 LMJP005 LMJF3 0 -9
JP LMJP005 0 0 1 -9
PJ LMPF2 0 0 2 -9
PJ LMJM2 0 0 1 -9
PJ LMPJ001 LMJM2 LMPF2 0 -9
PJ LMPJ002 LMJM2 LMPF2 0 -9
PJ LMPH003 LMJM2 LMPF2 0 -9
```

Then re-run the Mendelian analysis

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile all10_clean \
  --allow-extra-chr \
  --mendel \
  --out all10_clean_mendel
```

``` bash

cat all10_clean_mendel.imendel

 FID     IID   N
  JP LMJP005 42678
  JP   LMJF3 41898 
  JP   LMPM3 22438
  JP LMJP001 23576
  JP LMJP004 23288
  PJ   LMJM2 44286
  PJ   LMPF2 38861 
  PJ LMPJ001 23998
  PJ LMPJ002 23220
  PJ LMPH003 22445
```

The remaining Mendel errors are at the same scale in both families, so
they likely reflect genotype-call noise rather than pedigree/sample
identity problems.

Yes — that makes **substantial contamination less likely**.

If `LMPM3` were heavily contaminated or incorrectly merged, I’d expect
it to show elevated Mendel errors. Instead, after fixing the pedigree,
its count is right in the normal offspring range:

    LMPM3:   22,438
    LMJP001: 23,576
    LMJP004: 23,288

So the high sibling PI_HAT values are more likely due to
biology/estimation effects than a badly contaminated BAM.

I’d phrase the conclusion as:

> The corrected pedigree removes the major inconsistency, and Mendelian
> error rates are comparable across both families. This argues against a
> major BAM contamination or merge problem. The elevated PI_HAT among
> some JP siblings remains somewhat unusual, but it does not appear to
> be accompanied by excess Mendelian errors, so it is less concerning.

Still worth checking read groups, but this result is reassuring.

### Make bed file of regions with low mappability

I’m going to use the program GenMap to estimate mappability. I’ll try to
use the same filter Arne used for the bamsurgeon. Text from manuscript
draft: “First, we filtered the reference genome by estimating the
mappability across the genome using genmap v.X (REF) for 150bp kmers (-K
150) with up to 2 mismatches (-E 2) (150,2 mappability). We only kept
regions with unique mappability (value of 1). Second, we estimated the
coverage for uniquely mappable windows across the genome using mosdepth
v.X using the fast mode (–fast-mode) for each individual. We only used
windows with a coverage of at least 20x in the one assessed individual.”

I grabbed the genome from
`/fs/cbsubscb16/storage/silverside_progeny/genome/core/Mmenidia_refgenome_anchored.all_renamed_v2_core.fasta`.
Since I don’t have write access to that folder I copied it to my working
directory `/workdir/nina_july2026/silverside_progeny/genome`.

I then downloaded GenMap 1.3.0 with `conda install -c bioconda genmap`

Then I built the index

``` bash

cd /workdir/nina_july2026/silverside_progeny/genome

genmap index \
    -F Mmenidia_refgenome_anchored.all_renamed_v2_core.fasta \
    -I Mmenidia_refgenome_anchored.all_renamed_v2_core
```

Then compute mappability

``` bash

genmap map \
    -K 150 \
    -E 2 \
    -I Mmenidia_refgenome_anchored.all_renamed_v2_core \
    -O Mmenidia_refgenome_anchored.all_renamed_v2_core_mappability \
    -bg \
    -T 16
```

Convert the bedGraph to a BED of uniquely mappable regions (mappability
= 1) with:

``` bash

awk '$4==1 {print $1"\t"$2"\t"$3}' \
Mmenidia_refgenome_anchored.all_renamed_v2_core_mappability.bedgraph \
> Mmenidia_renamed_v2_core_uniquely_mappable.bed
```

Look at proportion retained

``` bash

awk '
FNR==NR {genome+=$2; next}
$4==1 {map+=$3-$2}
END {
    print "Uniquely mappable bases:", map
    print "Genome size:", genome
    print "Percent:", 100*map/genome
}' Mmenidia_refgenome_anchored.all_renamed_v2_core.fasta.fai \
Mmenidia_refgenome_anchored.all_renamed_v2_core_mappability.bedgraph
```

Uniquely mappable bases: 430587013

Genome size: 462853144

Percent: 93.0289

Next I’ll run mosdepth

``` bash

cd /workdir/nina_july2026/silverside_progeny/bams

mkdir -p mosdepth_1kb

for bam in *.bam; do
    sample=${bam%.bam}
    /programs/mosdepth-0.3.11/mosdepth -t 8 -b 1000 mosdepth_1kb/${sample}.1kb "$bam"
done
```

Then sum the read depth per 1kb window across all 10 individuals to look
at the distribution to determine cutoffs.

``` bash

cd /workdir/nina_july2026/silverside_progeny/bams/mosdepth_1kb

files=(*.regions.bed.gz)

zcat "${files[0]}" | cut -f1-3 > windows.tmp

for f in "${files[@]}"; do
    zcat "$f" | cut -f4 > "$(basename "$f").depth.tmp"
done

paste windows.tmp *.depth.tmp | \
awk '{
    sum=0
    for(i=4;i<=NF;i++) sum+=$i
    print $1"\t"$2"\t"$3"\t"sum
}' > all_bams.1kb.total_depth.bedgraph

rm windows.tmp *.depth.tmp
```

Then get quantiles

``` bash

cut -f4 all_bams.1kb.total_depth.bedgraph | sort -n | awk '
{a[NR]=$1}
END{
 print "N windows:", NR
 print "p1:", a[int(NR*0.01)]
 print "p5:", a[int(NR*0.05)]
 print "p50:", a[int(NR*0.50)]
 print "p95:", a[int(NR*0.95)]
 print "p99:", a[int(NR*0.99)]
 print "p99.5:", a[int(NR*0.995)]
}'
```

p1: 0 p5: 224.98 p50: 531.38 p95: 738.44 p99: 1747.28 p99.5: 2602.45

Based no this, I’m selected a min of 250 and max of 1750.

Then make a bad-depth BED

``` bash

LOW=250
HIGH=1750

awk -v low=$LOW -v high=$HIGH '$4 < low || $4 > high {print $1"\t"$2"\t"$3}' \
  all_bams.1kb.total_depth.bedgraph \
  > all_bams.1kb.bad_depth.bed
```

Interestingly, Áki had used a depth of 8000, but I think that might be
too high.

Check the distribution of combined read depth in the vcf

``` bash

bcftools query -f '%INFO/DP\n' all10.shared.clean.snps.vcf.gz > all10.shared.clean_site_total_depth.txt
```

Then summarize

``` bash

awk '
{sum+=$1; n++; if(min=="" || $1<min) min=$1; if($1>max) max=$1}
END {print "n="n, "mean="sum/n, "min="min, "max="max}
' all10.shared.clean_site_total_depth.txt
```

n=8314322 mean=1051.64 min=544 max=2492

The mean of ~1000 is much higher than the p50 of mosdepth, which was
~531. I wonder if the variant set is biased from too aggessive
filtering, or many lower-covered windows in the mosdepth data got
filtered out.

## Filter out SNPs by mappability and mosdepth

``` bash

echo "Starting SNPs:"
bcftools view -H all10.shared.clean.snps.vcf.gz | wc -l
8314322

echo "After both mappability + mosdepth:"
bcftools view \
  -R ../genome/Mmenidia_renamed_v2_core_uniquely_mappable.bed \
  -T ^../bams/mosdepth_1kb/all_bams.1kb.bad_depth.bed \
  -H all10.shared.clean.snps.vcf.gz | wc -l
19642
```

### Filter out inversions

There’s a list of inversions Azwad has used here:
`/workdir/azwad/haplotagging_ms/docs/silverside_all_sv_regions_sorted.bed`

This BED file uses `chr01`, but our VCF earlier used names like
`Mme_chr06`. The chromosome names must match exactly, so I’ll add the
Mme\_ prefix to all the chrom names

``` bash

awk 'BEGIN{OFS="\t"} NR==1 {print; next} {$1="Mme_"$1; print}' \
/workdir/azwad/haplotagging_ms/docs/silverside_all_sv_regions_sorted.bed \
> sv_regions_Mme.bed
```

## Filter the vcf based on mappability \< 1, exclude inversions, and filter by mosdepth.

First merge the exclusion beds

``` bash

# There were issues with the bed file formatting, so re-writing both

awk 'BEGIN{OFS="\t"} 
     NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' \
../bams/mosdepth_1kb/all_bams.1kb.bad_depth.bed \
> ../bams/mosdepth_1kb/all_bams.1kb.bad_depth.clean.bed

awk 'BEGIN{OFS="\t"} 
     NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1,$2,$3}' \
../genome/sv_regions_Mme.noheader.bed \
> ../genome/sv_regions_Mme.clean.bed

cat ../bams/mosdepth_1kb/all_bams.1kb.bad_depth.clean.bed \
    ../genome/sv_regions_Mme.clean.bed | \
sort -k1,1 -k2,2n | \
bedtools merge > ../genome/mosdepth_and_sv_exclude_regions.bed
```

Now filter. I’m going to also remove sites for which any individual has
a depth \<20. The rreason is that you’re interested in **high-confidence
Mendelian inheritance and de novo mutation detection**. A missing or
low-depth genotype in even one individual (especially a parent) makes
the site much less informative. Rather than setting that genotype to
missing and carrying around partially observed sites, it’s cleaner to
remove the site altogether.

``` bash

bcftools view \
    -R ../genome/Mmenidia_renamed_v2_core_uniquely_mappable.bed \
    -T ^../genome/mosdepth_and_sv_exclude_regions.bed \
    -m2 -M2 \
    -v snps \
    -f PASS \
    -i 'MIN(FMT/DP)>=20' \
    all10.shared.clean.snps.vcf.gz \
    -Oz \
    -o all10.shared.clean.mappability_depth_filtered.vcf.gz
    
    
# Count retained SNPs
bcftools view -H all10.shared.clean.mappability_depth_filtered.vcf.gz | wc -l
6441729



tabix -p vcf all10.shared.clean.mappability_depth_filtered.vcf.gz
    
```

## Re-do the relatedness and Mendelian error analysis with the filtered dataset

## Convert to PLINK

I want to add SNP positions as the name for each site so I can track
where

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --vcf all10.shared.clean.mappability_depth_filtered.vcf.gz \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --make-bed \
  --out all10.shared.clean.mappability_depth_filtered
```

## Edit the `.fam` file (created by plink in the conversion)

I use the updated pedigree

FamilyID IndividualID FatherID MotherID Sex Phenotype

``` bash
nano all10.shared.clean.mappability_depth_filtered.fam

JP LMJF3 0 0 2 -9
JP LMPM3 LMJP005 LMJF3 0 -9
JP LMJP001 LMJP005 LMJF3 0 -9
JP LMJP004 LMJP005 LMJF3 0 -9
JP LMJP005 0 0 1 -9
PJ LMPF2 0 0 2 -9
PJ LMJM2 0 0 1 -9
PJ LMPJ001 LMJM2 LMPF2 0 -9
PJ LMPJ002 LMJM2 LMPF2 0 -9
PJ LMPH003 LMJM2 LMPF2 0 -9
```

## Run Mendelian error analysis

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile all10.shared.clean.mappability_depth_filtered \
  --allow-extra-chr \
  --mendel \
  --out all10.shared.clean.mappability_depth_filtered_mendel
```

Useful outputs

``` bash

cat all10.shared.clean.mappability_depth_filtered_mendel.imendel

FID     IID   N
  JP LMJP005 31379
  JP   LMJF3 33866 
  JP   LMPM3 17353
  JP LMJP001 17877
  JP LMJP004 17770
  PJ   LMJM2 32613
  PJ   LMPF2 31404 
  PJ LMPJ001 18365
  PJ LMPJ002 17862
  PJ LMPH003 17062
```

## Estimate relatedness among all 10 individuals

``` bash

/programs/plink-1.9-x86_64-beta7/plink \
  --bfile all10.shared.clean.mappability_depth_filtered \
  --allow-extra-chr \
  --genome \
  --out all10.shared.clean.mappability_depth_filtered_relatedness
```

``` bash

for chr in $(seq -w 1 24); do
  c="Mme_chr${chr}"

  bcftools view -r $c all10.shared.clean.mappability_depth_filtered.vcf.gz \
    -Oz -o ${c}.vcf.gz
  tabix -p vcf ${c}.vcf.gz

  n=$(bcftools view -H ${c}.vcf.gz | wc -l)

  /programs/plink-1.9-x86_64-beta7/plink \
    --vcf ${c}.vcf.gz \
    --allow-extra-chr \
    --make-bed \
    --out ${c}

  /programs/plink-1.9-x86_64-beta7/plink \
    --bfile ${c} \
    --allow-extra-chr \
    --genome \
    --out ${c}

  awk -v chr=$c -v n=$n '
    NR>1 {
      pair=$2"-"$4
      print chr, n, $2, $4, $10
    }' ${c}.genome
done > per_chrom_PIHAT.tsv


awk '$3=="LMPM3" && ($4=="LMJP001" || $4=="LMJP004") || \
     $4=="LMPM3" && ($3=="LMJP001" || $3=="LMJP004") || \
     $3=="LMJP001" && $4=="LMJP004" || \
     $3=="LMJP004" && $4=="LMJP001"' \
per_chrom_PIHAT.tsv
```

This didn’t get to the bottom of things. Will park this for now.

## Discover de novo SNPs from the family-specific vcf

### For one offspring: `LMPJ001`

This finds sites where both parents are the same homozygote, `LMPJ001`
is heterozygous, and siblings are not heterozygous.

I’ll use the filters Áki describes in the manuscript

\`Only variants with an allelic depth (AD) ratio between 0.25 and .75
were considered. Also, candidate de novo variants in a progeny that had
any level of minor allele AD value in either parent or a sibling were
not considered. Only heterozygous variants observed in a single
offspring, and not in either parent or sibling, were considered as
candidate de novo mutations.

BCFtools results had the following filters applied: Mann-Whitney
rank-sum test of Mapping Quality Bias (MQB) \< 0.06, Mann-Whitney
rank-sum test of Read Position Bias (RPB) \< 0.06, Mann-Whitney rank-sum
test of Base Quality Bias (BQB) \< 0.06 (MQB, RPB, and BQB are unique
fields to BCFtools). Both BCFtools and GATK4 results were filtered with
Quality by Depth (QD)\<2.0, and average Mapping Quality (MQ) \<40.0.
Additionally the following filters were applied to GATK4 results: strand
bias (FS)\>60.0, Strand Odds Ratio (SOR)\>3.0,  Mann-Whitney rank-sum
test of Mapping Quality (MQRankSum) \< -12.5, Mann-Whitney rank-sum test
of site position within reads (ReadPosRankSum) \< -8.0. \`

Then for consistency with bamsurgeon,

``` bash

VCF=PJ_bwa_GATK4_miss_g5_Fullfilt_FRMTfilt.snps.vcf.gz

for target in LMPJ001 LMPJ002 LMPH003; do

  if [ "$target" = "LMPJ001" ]; then sib1=LMPJ002; sib2=LMPH003; fi
  if [ "$target" = "LMPJ002" ]; then sib1=LMPJ001; sib2=LMPH003; fi
  if [ "$target" = "LMPH003" ]; then sib1=LMPJ001; sib2=LMPJ002; fi

  bcftools view \
    -m2 -M2 -v snps \
    -i 'QD>=2.0 && MQ>=40.0 && FS<=60.0 && SOR<=3.0 && MQRankSum>=-12.5 && ReadPosRankSum>=-8.0' \
    "$VCF" | \
  bcftools query \
    -s LMPF2,LMJM2,$target,$sib1,$sib2 \
    -f '%CHROM\t%POS[\t%GT\t%DP\t%AD]\n' | \
  awk 'BEGIN{OFS="\t"}
  function hom(g){return g=="0/0" || g=="0|0" || g=="1/1" || g=="1|1"}
  function het(g){return g=="0/1" || g=="1/0" || g=="0|1" || g=="1|0"}
  function refAD(ad){split(ad,a,","); return a[1]+0}
  function altAD(ad){split(ad,a,","); return a[2]+0}
  function ab_ref(ad){split(ad,a,","); return (a[1]+a[2]>0) ? a[1]/(a[1]+a[2]) : 0}
  function ab_alt(ad){split(ad,a,","); return (a[1]+a[2]>0) ? a[2]/(a[1]+a[2]) : 0}
  $3==$6 && hom($3) && het($9) && !het($12) && !het($15) && $4>=20 && $7>=20 && $10>=20 && $13>=20 && $16>=20 && ((($3=="0/0" || $3=="0|0") && ab_alt($11)>=0.25 && ab_alt($11)<=0.75 && altAD($5)==0 && altAD($8)==0 && altAD($14)==0 && altAD($17)==0) || (($3=="1/1" || $3=="1|1") && ab_ref($11)>=0.25 && ab_ref($11)<=0.75 && refAD($5)==0 && refAD($8)==0 && refAD($14)==0 && refAD($17)==0)) {print $1, $2-1, $2}
  ' > ${target}_cDeNovo_candidates.bed

  if [ -s ${target}_cDeNovo_candidates.bed ]; then
    bcftools view \
      -R ${target}_cDeNovo_candidates.bed \
      "$VCF" \
      -Oz \
      -o ${target}_cDeNovo_candidates.vcf.gz

    tabix -p vcf ${target}_cDeNovo_candidates.vcf.gz
  fi

  echo "$target:"
  wc -l ${target}_cDeNovo_candidates.bed

done
```

## Remaining things to do

- Figure out if we need to account for the higher relatedness among some
  inds

- Assess whether the Mendelian inconsistencies occur in certain genomic
  regions we can filter out. The list of Mendelian inconsistencies are
  in
  `/workdir/nina_july2026/silverside_progeny/vcf/all10_clean_mendel.mendel`

- Figure out where the sample swap occurred, and if it was
  bioinformatically, how far we need to go back to correct it

- Figure out why Áki’s candidate de novos don’t match the manual
  inspection, even for PJ where there are no pedigree issues

## Checking read groups (to look for evidence of reads from different individuals getting merged)

Nothing suspicious detected here, but we’ll have to look at the approach
used for mergining to see if read group names were added manually so
this could mask true mix-up of bam files from different individuals

``` bash

samtools view -H LMPM3_bwa.sorted.RG.MD.merged.bam | grep '^@RG'

@RG ID:HC3LCDSX2.L1 LB:Baym_HC3LCDSX2.L1    PL:ILLUMINA SM:LMPM3    PU:HC3LCDSX2.L1.A00261.402
@RG ID:HC553DSX2.L1 LB:Baym_HC553DSX2.L1    PL:ILLUMINA SM:LMPM3    PU:HC553DSX2.L1.A00261.408
@RG ID:HVCVFDSXY.L3 LB:Baym_HVCVFDSXY.L3    PL:ILLUMINA SM:LMPM3    PU:HVCVFDSXY.L3.A00261.357
```

``` bash
samtools view -H LMJP001_bwa.sorted.RG.MD.merged.bam | grep '^@RG'

@RG ID:HC3LCDSX2.L1 LB:Baym_HC3LCDSX2.L1    PL:ILLUMINA SM:LMJP001  PU:HC3LCDSX2.L1.A00261.402
@RG ID:HC553DSX2.L1 LB:Baym_HC553DSX2.L1    PL:ILLUMINA SM:LMJP001  PU:HC553DSX2.L1.A00261.408
@RG ID:HVCVFDSXY.L3 LB:Baym_HVCVFDSXY.L3    PL:ILLUMINA SM:LMJP001  PU:HVCVFDSXY.L3.A00261.357
```

``` bash

samtools view -H LMJP004_bwa.sorted.RG.MD.merged.bam | grep '^@RG'

@RG ID:HC3LCDSX2.L1 LB:Baym_HC3LCDSX2.L1    PL:ILLUMINA SM:LMJP004  PU:HC3LCDSX2.L1.A00261.402
@RG ID:HC553DSX2.L1 LB:Baym_HC553DSX2.L1    PL:ILLUMINA SM:LMJP004  PU:HC553DSX2.L1.A00261.408
@RG ID:HVCVFDSXY.L3 LB:Baym_HVCVFDSXY.L3    PL:ILLUMINA SM:LMJP004  PU:HVCVFDSXY.L3.A00261.357
```
