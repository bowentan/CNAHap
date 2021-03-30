# CNAHap
CNAHap is a tool to estimate tumor allele-specific copy numbers and perform germline
haplotyping along segments with imbalanced copy numbers.

## Download and installation
Please make sure you have Python3.6 or above to succussfully run CNAHap. 

To get CNAHap, you can clone this repository to your workspace
```bash
git clone git@github.com:deepomicslab/CNAHap.git
```

## Usage
```
usage: cnahap.py [-h] -f CNV_FILE -n NORMAL_VCF -t TUMOR_VCF [-v MIN_NUM_HET_VARIANTS] -p PURITY -o OUT_FN -s SAMPLENAME -c OUT_CHECK

optional arguments:
  -h, --help            show this help message and exit
  -f CNV_FILE, --cnv-file CNV_FILE
                        CNV segment file from Accurity
  -n NORMAL_VCF, --normal-vcf NORMAL_VCF
                        Path of called VCF for normal cells
  -t TUMOR_VCF, --tumor-vcf TUMOR_VCF
                        Path of called VCF for tumor cells
  -v MIN_NUM_HET_VARIANTS, --min-num-het-variants MIN_NUM_HET_VARIANTS
                        Minimum number of heterozygeous SNVs in a segment. Default: 30
  -p PURITY, --purity PURITY
                        Tumor purity. Default: 1
  -o OUT_FN, --out-vcf OUT_FN
                        Phased output file path in VCF
  -s SAMPLENAME, --sample SAMPLENAME
                        Sample name
  -c OUT_CHECK, --check-file OUT_CHECK
                        Check results/Copy number estimate results
```

## Output for copy number estimation (-c)
One example is provided below.
```
chrom start end accurity_cp accurity_major_cp n_snps binom_ais avg_depths M_depths m_depths ad_diff dp_phased cp M m phased
chr1 1 13000 na na 0 na na na na na False na na na False
chr1 13001 6840000 2.0 1.0 4431 0.29 32.76 18.71 14.05 4.66 False na na na False
chr1 6840001 7834500 na na 954 0.67 129.56 96.67 32.89 63.77 True 8.0 6.0 2.0 True
chr1 7834501 18952000 2.0 1.0 7072 0.28 32.68 18.63 14.05 4.58 False na na na False
chr1 18952001 18952500 na na 0 na na na na na False na na na False
chr1 18952501 19080000 4.0 2.0 83 0.18 66.75 36.43 30.32 6.10 True 4.0 2.0 2.0 False
chr1 19080001 19080500 na na 1 0.36 77.0 47.0 30.0 17.0 False na na na False
chr1 19080501 121099000 2.0 1.0 58900 0.29 32.62 18.62 14.00 4.61 False na na na False
chr1 121099001 121132000 na na 4 0.52 29.25 18.75 10.5 8.25 False na na na False
```

The explanations for the metrics in the header is:

- chrom: the chrom/contig of the segment
- start: the start position of the segment
- end: the end position of the segment
- accurity_cp: the total copy numbers from Accurity
- accurity_major_cp: the major copy numbers
- n_snps: the number of SNPs in the segment
- binom_ais: the allele imbalance value normalized by binomial coefficients
- avg_depths: the average depth of SNPs
- M_depths: the average depth of alleles with the larger depths
- m_depths: the average depth of alleles with the smaller depths
- ad_diff: the average allele difference
- dp_phased: whether the segment is considered as normall segment with 2 copies
- cp: the total copy numbers from CNAHap
- M: the major copy from CNAHap
- m: the minor copy from CNAHap
- phased: whether the segment is considered to be phased

## Citation
Please cite our paper if you use CNAHap in your work.

Bowen Tan+, Lingxi Chen+, Wenlong Jia, Yanfei Wang, Hechen Li, Shuai Cheng Li*, CNAHap: a germline haplotyping method using tumor allele-specific copy number alteration, bioRxiv 2021.03.27.437314; doi: https://doi.org/10.1101/2021.03.27.437314