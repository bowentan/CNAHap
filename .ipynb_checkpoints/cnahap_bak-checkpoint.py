import sys
import argparse
import gzip
import numpy as np
import scipy as sp
from pysam import VariantFile
import pandas as pd
import subprocess
import sklearn.cluster as cl
from Bio import bgzf
import mip

# def get_region_recs(vcf, region, sample):
# #     print(region)
#     recs = np.array([
#         (r.pos, r.ref, r.alts[0], *r.samples[sample]['GT'], *r.samples[sample]['AD'])
#         for r in vcf.fetch(region[0], int(region[1]) - 1, int(region[2]))
#         if len(r.ref) == 1 and len(r.alts) == 1 and len(r.alts[0]) == 1
#     ])
# #     print(recs)
#     return recs

def segmentation(regs):
    res = []
    for c in [f'chr{i}' for i in range(1, 23)]:
        c_regs = regs[regs['chrom'] == c]
        bps = np.sort(np.hstack([c_regs['start'], c_regs['end']]))
        
        start = 1
        for p in bps:
            res.append((c, start, p - 1))
            start = p
        res.append((c, start, -1))
    return res
        

def get_overlap(reg1, reg2):
    if reg1[0] != reg2[0]:
        return None
    start = max(reg1[1], reg2[1])
    end = min(reg1[2], reg2[2])
    if start > end:
        return None
    return reg1[0], start, end
    
def get_recs(vcf_name, regs):
    cmd = [
        '/home/BIOINFO_TOOLS/mutation_tools/bcftools/bcftools-1.9/bcftools', 
        'view', 
        vcf_name, 
        '-H',
        '-g', 
        'het',
        '-v',
        'snps',
        '-m2',
        '-M2',
        '-r',
        ','.join([r[2] > 0 and f'{r[0]}:{r[1]}-{r[2]}' or f'{r[0]}:{r[1]}-' for r in regs])
    ]
    recs = []
    for r in subprocess.check_output(cmd).decode('utf-8')[:-1].split('\n'):
        r = r.split('\t')
        ad = list(map(int, r[-1].split(':')[1].split(',')))
        recs.append((r[0], int(r[1]), r[3], r[4], *ad))
    dtypes = [('chrom', '<U5'), ('pos', int), 
              ('ref', '<U1'), ('alt', '<U1'), 
              ('ad1', int), ('ad2', int)]
    return np.array(recs, dtype=dtypes)
    
# def get_region_recs(vcf, region, sample):
# #     print(region)
#     recs = np.array([
#         (r.pos, r.ref, r.alts[0], *r.samples[sample]['GT'], *r.samples[sample]['AD'])
#         for r in vcf.fetch(region[0], region[1] - 1, region[2])
#         if len(r.ref) == 1 and len(r.alts) == 1 and len(r.alts[0]) == 1
#     ], dtype=[('pos', int), ('ref', '<U1'), ('alt', '<U1'), ('gt1', int), ('gt2', int), ('ad1', int), ('ad2', int)])
# #     print(recs)
#     return recs

def get_region_recs(recs, region):
    return recs[(recs['chrom'] == region[0]) & (recs['pos'] >= region[1]) & (recs['pos'] <= region[2])]
    
def filter_recs(normal_recs, tumor_recs):
    # Drop non-het
#     normal_recs = normal_recs[normal_recs['gt1'] + normal_recs['gt2'] == 1]
#     tumor_recs = tumor_recs[tumor_recs['gt1'] + tumor_recs['gt2'] == 1]
    
    # Intersection between normal and tumor
    fields = ['chrom', 'pos', 'ref', 'alt']
    inter, normal_inter_idx, tumor_inter_idx = np.intersect1d(normal_recs[fields], tumor_recs[fields], return_indices=True)
    normal_recs = normal_recs[normal_inter_idx]
    tumor_recs = tumor_recs[tumor_inter_idx]
    
    # Compare REF and ALT
#     is_equal = normal_recs[['ref', 'alt']] == tumor_recs[['ref', 'alt']]
#     normal_recs = normal_recs[is_equal]
#     tumor_recs = tumor_recs[is_equal]
    
    # Drop ad == 0
    tumor_recs = tumor_recs[np.logical_and(tumor_recs['ad1'] > 0, tumor_recs['ad2'] > 0)]
    return tumor_recs
    
def get_ad_old(normal_vcf, tumor_vcf, regions, sample_normal, sample_tumor):
    orig_ads = []
    for region in regions:
        orig_ad_l = []
        de_ad_l = []
        for normal_record in normal_vcf.fetch(region[0], int(region[1]) - 1, int(region[2])):
            if sum(normal_record.samples[sample_normal]['GT']) != 1:
    #             print(f"skip GT, {normal_record.chrom} {normal_record.pos} {normal_record.samples[sample_normal]['GT']}")
                continue
            if len(normal_record.ref) > 1 or any(map(lambda x : len(x) > 1, normal_record.alts)):
                # Skip the non-SNV loci
    #             print(f'skip non-SNV, {normal_record.CHROM} {normal_record.POS} {normal_record.ref} {normal_record.alts}')
                continue
            try:
                tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
                if normal_record.ref != tumor_record.ref or normal_record.alts != tumor_record.alts:
                    continue
                if sum(tumor_record.samples[sample_tumor]['GT']) != 1:
                    continue
                ad = np.array(tumor_record.samples[sample_tumor]['AD'])
                if any(ad <= 0):
                    continue
#                 de_ad = ad - round(sum(ad) * (1 - frac) / 2)
#                 if any(de_ad <= 0):
#                     continue
                orig_ad_l.append(ad)
#                 de_ad_l.append(de_ad)
            except:
                pass
        orig_ad_l = np.array(orig_ad_l)
#         de_ad_l = np.array(de_ad_l)
#         print(f'{region[0]} {region[1]} {region[2]} {orig_ad_l.size}', file=sys.stderr)
        orig_ads.append(orig_ad_l)
#         de_ads.append(de_ad_l)
#         n_snps.append(orig_ad_l.size)
    return np.array(orig_ads)

def get_region_ads(tumor_recs):
    return np.array(tumor_recs[['ad1', 'ad2']].tolist())

def cal_region_ads_by_frac(recs, major_frac, minor_frac):
    if recs.size == 0:
        return np.empty(0, dtype=int)
    ads = get_region_ads(recs)
    ads_major = ads.max(axis=1)
    ads_minor = ads.min(axis=1)
    res_ads = np.array(list(zip((ads_major * major_frac).round(), (ads_minor * minor_frac).round())))
    return res_ads[(res_ads > 0).sum(axis=1) == 2]

def cal_tumor_dna_frac(regions, purity):
    M, m = regions['M'], regions['m']
    tumor_major_frac = M * purity / ((1 - purity) * 1 + purity * M)
    tumor_minor_frac = m * purity / ((1 - purity) * 1 + purity * m)
    if isinstance(tumor_major_frac, pd.Series):
        return tumor_major_frac.values, tumor_minor_frac.values
    return tumor_major_frac, tumor_minor_frac

def cal_region_avg_depth(ads):
    if ads.size == 0:
        return np.nan
    return ads.sum(axis=1).mean()

def cal_region_avg_M(ads):
    if ads.size == 0:
        return np.nan
    return ads.max(axis=1).mean()

def cal_region_avg_m(ads):
    if ads.size == 0:
        return np.nan
    return ads.min(axis=1).mean()

def cal_region_ad_diff(ads):
    if ads.size > 0:
        return np.absolute(np.diff(ads, axis=1)).mean()
    else:
        return np.nan

def cal_region_ad_diff_std(ads):
    if ads.size > 0:
        return np.absolute(np.diff(ads, axis=1)).std()
    else:
        return np.nan

def cal_log10_binom_pmf(k, n, p):
    coef_nom_log10 = np.log10(np.arange(1, n + 1)).sum()

    coef_denom_log10 = np.log10(np.arange(1, k + 1)).sum() + np.log10(np.arange(1, n - k + 1)).sum()
    return coef_nom_log10 - coef_denom_log10 + k * np.log10(p) + (n - k) * np.log10(1 - p)

def cal_region_binom_ais_mean(ads):
    if ads.size == 0:
        return np.nan
    vfunc = np.vectorize(cal_log10_binom_pmf)
    M = ads.max(axis=1)
    m = ads.min(axis=1)
    w = (np.log10(2) + vfunc(M, M + m, 0.5)) * -10
    return (w * (M - m) / M).sum() / w.sum()
#     if ad_l.size == 0:
#         return np.nan
#     w = []
#     s = []
#     for ad in ad_l:
#         M = max(ad)
#         m = min(ad)
# # 		cw = np.log10(sp.stats.binom.pmf(M, M + m, 0.5) * 2) * (-10)
#         cw = (np.log10(2) + cal_log10_binom_pmf(M, M + m, 0.5)) * (-10)
#         w.append(cw)
#         if np.isnan(cw * (M - m) / M):
#             print(f'{M} {m}')
#         s.append(cw * (M - m) / M)
#     return sum(s) / sum(w)

def check_normal(df, min_n_snp):
    # depth cluster
    cp = np.full(len(df), 0)
#     df_non0 = df.loc[lambda row: row.n_snps > 0].loc[lambda row: row.m > 0]
    df_non0 = df.loc[lambda row: row.n_snps > min_n_snp]
    tot_avg_depth = (df_non0.avg_depths * df_non0.n_snps).sum() / df_non0.n_snps.sum() / 2
    dp2cp = df_non0.avg_depths / tot_avg_depth
    cp[df_non0.index] = dp2cp
    # n_cluster = int(np.ceil(dp2cp.max() - dp2cp.min())) - 1
    init_centers = np.arange(round(dp2cp.min()), round(dp2cp.max()) + 1).reshape(-1, 1) * tot_avg_depth
    n_cluster = init_centers.size
    # n_cluster = int(round(dp2cp.max())) - 1
    dp_ms = cl.KMeans(n_cluster, init_centers, n_init=1, max_iter=1).fit(df_non0.avg_depths.values.reshape(-1, 1))
        
    df_dp = df_non0.loc[dp_ms.labels_ != dp_ms.cluster_centers_.argmin()]
    dp_phased = np.full(len(df), False)
    dp_phased[df_dp.index] = True
    
    res = pd.concat([df, pd.DataFrame({'dp_phased' : dp_phased})], axis=1)
    return dp_ms, res

def check_balance(df, min_n_snp):
    # depth cluster
#     cp = np.full(len(df), 0)
# #     df_non0 = df.loc[lambda row: row.n_snps > 0].loc[lambda row: row.m > 0]
#     df_non0 = df.loc[lambda row: row.n_snps > min_n_snp]
#     tot_avg_depth = (df_non0.avg_depths * df_non0.n_snps).sum() / df_non0.n_snps.sum() / 2
#     dp2cp = df_non0.avg_depths / tot_avg_depth
#     cp[df_non0.index] = dp2cp
#     # n_cluster = int(np.ceil(dp2cp.max() - dp2cp.min())) - 1
#     init_centers = np.arange(round(dp2cp.min()), round(dp2cp.max()) + 1).reshape(-1, 1) * tot_avg_depth
#     n_cluster = init_centers.size
#     # n_cluster = int(round(dp2cp.max())) - 1
#     dp_ms = cl.KMeans(n_cluster, init_centers, n_init=1, max_iter=1).fit(df_non0.avg_depths.values.reshape(-1, 1))
#     # min_cluster_index = df_non0.loc[dp_ms.labels_ == dp_ms.cluster_centers_.argmin()].index
#     # orig_index = np.full(len(min_cluster_index), -1)
#     # while orig_index.size != min_cluster_index.size or not all(min_cluster_index == orig_index):
#     #     orig_index = min_cluster_index
        
#     #     normal_dp_like = df_non0.loc[dp_ms.labels_ == dp_ms.cluster_centers_.argmin()]
#     #     tot_avg_depth = (normal_dp_like.avg_depths * normal_dp_like.n_snps).sum() / normal_dp_like.n_snps.sum() / 2
#     #     dp2cp = df_non0.avg_depths / tot_avg_depth
#     # #     n_cluster = int(np.ceil(dp2cp.max() - dp2cp.min())) - 1
#     #     init_centers = np.arange(dp2cp.min(), dp2cp.max()).reshape(-1, 1)
#     #     n_cluster = init_centers.size
#     # #     n_cluster = int(round(dp2cp.max())) - 1
#     #     dp_ms = cl.KMeans(n_cluster, max_iter=1, n_init=1).fit(dp2cp.values.reshape(-1, 1))
#     #     min_cluster_index = df_non0.loc[dp_ms.labels_ == dp_ms.cluster_centers_.argmin()].index
        
#     df_dp = df_non0.loc[dp_ms.labels_ != dp_ms.cluster_centers_.argmin()]
#     dp_phased = np.full(len(df), False)
#     dp_phased[df_dp.index] = True
    
    # ais cluster
    df_dp = df.loc[lambda row: row.dp_phased == True]
    df_n = df_dp.loc[lambda row: row.n_snps > min_n_snp]
    aisXstds_ms = cl.KMeans(n_cluster).fit(df_n.aisXstds.values.reshape(-1, 1))
    ais_ms = cl.KMeans(n_cluster).fit(df_n.binom_ais.values.reshape(-1, 1))
    # res = df_n.loc[df_n.index.isin(df_dp.loc[ais_ms.labels_ != ais_ms.cluster_centers_.argmin()].index)]
    res = df_n.loc[np.logical_and(ais_ms.labels_ != ais_ms.cluster_centers_.argmin(), aisXstds_ms.labels_ != aisXstds_ms.cluster_centers_.argmin())]
    phased = np.full(len(df), False)
    phased[res.index] = True
#     res = pd.concat([df, pd.DataFrame({'cp' : cp, 'dp_phased' : dp_phased, 'phased' : phased})], axis=1)
    res = pd.concat([df, pd.DataFrame({'phased' : phased})], axis=1)
    return ais_ms, aisXstds_ms, res

def estimate_cp(bal_res, de_ads):
    df_normal = bal_res.loc[lambda row: row.dp_phased == False]
    df_abnormal = bal_res.loc[lambda row: row.dp_phased == True]
    ads_abnormal = de_ads[df_abnormal.index]
    
    M_, m_ = np.full(len(bal_res), np.nan), np.full(len(bal_res), np.nan)
    hap_dp = df_normal.avg_depths.mean() / 2
    for row in df_abnormal.itertuples():
        ads = de_ads[row.Index]
        mod = mip.Model()
        M, m = mod.add_var(var_type=mip.INTEGER, lb=1), mod.add_var(var_type=mip.INTEGER, lb=1)
        eM = [mod.add_var() for i in range(len(ads))]
        em = [mod.add_var() for i in range(len(ads))]

        for i in range(len(ads)):
            mod += ads[i].max() - M * hap_dp >= -eM[i]
            mod += ads[i].max() - M * hap_dp <= eM[i]

            mod += ads[i].min() - m * hap_dp >= -em[i]
            mod += ads[i].min() - m * hap_dp <= em[i]
        mod += M >= m
        mod.objective = mip.xsum(eM[i] for i in range(len(ads))) + mip.xsum(em[i] for i in range(len(ads)))
        mod.optimize()
        M_[row.Index] = M.x
        m_[row.Index] = m.x
    res = pd.concat([
        bal_res,
        pd.DataFrame({
            'cp': M_ + m_,
            'M': M_,
            'm': m_
        })
    ], axis=1)
    return res

# def phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, purity, sample_normal, sample_tumor):
#     segment_first = True
#     num = 0
#     for normal_record in normal_vcf.fetch(region[0], int(region[1]) - 1, int(region[2])):
#         if sum(normal_record.samples[sample_normal]['GT']) != 1:
#             # filter out homozygotes
#             continue
#         if len(normal_record.ref) > 1 or any(map(lambda x : len(x) > 1, normal_record.alts)):
#             # Skip the non-SNV loci
#             continue
#         # catch exception that the record doesn't exist
#         try:
#             tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
#             if normal_record.ref != tumor_record.ref or normal_record.alts != tumor_record.alts:
#                 continue
#             if sum(tumor_record.samples[sample_tumor]['GT']) != 1:
#                 continue
#             ad = list(map(lambda x: round(x * purity), tumor_record.samples[sample_tumor]['AD']))
#             if any(ad) == 0:
#                 continue
#             # if sum(ad) < 5:
#             #     continue
#             if ad[0] > ad[1]:
#                 tumor_record.samples[sample_tumor]['GT'] = (0, 1)
#                 tumor_record.samples[sample_tumor]['PB'] = np.log10(sp.stats.binom.pmf(ad[0], sum(ad), 0.5) * 2) * (-10)
#             else:
#                 tumor_record.samples[sample_tumor]['GT'] = (1, 0)
#                 tumor_record.samples[sample_tumor]['PB'] = np.log10(sp.stats.binom.pmf(ad[1], sum(ad), 0.5) * 2) * (-10)

#             if segment_first:
#                 tumor_record.samples[sample_tumor].phased = False
#             else:
#                 tumor_record.samples[sample_tumor].phased = True
#             phased_vcf.write(tumor_record)
#             num += 1
#             segment_first = False
#         except:
#             pass

# def phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, purity, sample_normal, sample_tumor):
#     segment_first = True
#     num = 0
#     for normal_record in normal_vcf.fetch(region[0], int(region[1]) - 1, int(region[2])):
#         if sum(normal_record.samples[sample_normal]['GT']) != 1:
#             # filter out homozygotes
#             continue
#         if len(normal_record.ref) > 1 or any(map(lambda x : len(x) > 1, normal_record.alts)):
#             # Skip the non-SNV loci
#             continue
#         # catch exception that the record doesn't exist
#         try:
#             tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
#             if normal_record.ref != tumor_record.ref or normal_record.alts != tumor_record.alts:
#                 continue
#             if sum(tumor_record.samples[sample_tumor]['GT']) != 1:
#                 continue
#             ad = list(map(lambda x: round(x * purity), tumor_record.samples[sample_tumor]['AD']))
#             if any(ad) == 0:
#                 continue
#             # if sum(ad) < 5:
#             #     continue
#             if ad[0] > ad[1]:
#                 tumor_record.samples[sample_tumor]['GT'] = (0, 1)
#                 tumor_record.samples[sample_tumor]['PB'] = np.log10(sp.stats.binom.pmf(ad[0], sum(ad), 0.5) * 2) * (-10)
#             else:
#                 tumor_record.samples[sample_tumor]['GT'] = (1, 0)
#                 tumor_record.samples[sample_tumor]['PB'] = np.log10(sp.stats.binom.pmf(ad[1], sum(ad), 0.5) * 2) * (-10)

#             if segment_first:
#                 tumor_record.samples[sample_tumor].phased = False
#             else:
#                 tumor_record.samples[sample_tumor].phased = True
#             phased_vcf.write(tumor_record)
#             num += 1
#             segment_first = False
#         except:
#             pass
        
def phase_segment(out_filename, res, tumor_recs_used_reg, de_ads):
    with bgzf.BgzfWriter(out_filename) as fout:
        for info, ads in zip(tumor_recs_used_reg[res.phased], de_ads[res.phased]):
            segment_first = True
            first_max = -1
            for rec, ad in zip(info, ads):
                line = f'{rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}\t{rec[4]}\t{rec[5]}\t{ad[0]}\t{ad[1]}'
                if first_max < 0:
                    if ad[0] > ad[1]:
                        first_max = 0
                    else:
                        first_max = 1
                if segment_first:
                    if first_max == 0:
                        if ad[0] > ad[1]:
                            line += '\t0/1'
                        else:
                            line += '\t1/0'
                    else:
                        if ad[0] > ad[1]:
                            line += '\t1/0'
                        else:
                            line += '\t0/1'
                    segment_first = False
                else:
                    if first_max == 0:
                        if ad[0] > ad[1]:
                            line += '\t0|1'
                        else:
                            line += '\t1|0'
                    else:
                        if ad[0] > ad[1]:
                            line += '\t1|0'
                        else:
                            line += '\t0|1'
                fout.write(line + '\n')
            
def region_parser(args):
    if args.bed_file is not None:
        # BED file is given
        regions = []
        try:
            fin = gzip.open(args.bed_file, 'rb')
            line = fin.readline()[:-1]
        except:
            fin = open(args.bed_file, 'r')
            line = fin.readline()[:-1]
        while line:
            line = line.split()
            region = line[0:3]
            chrom = region[0]
            start, end = map(int, region[1:3])
            total_copy, minor = map(int, line[3].split('m'))
            major = total_copy - minor
            if start > end:
                raise ValueError(f'invalid region: chrom ({chrom}), ({start}) > end ({end})')
            regions.append((chrom, start, end, major, minor))
            # regions.append([chrom, start, end])
            line = fin.readline()[:-1]
#     else:
#         chrom = args.region[0].split(':')[0]
#         start, end = map(int, args.region.split(':')[1].split('-'))
#         total_copy, minor = map(int, args.region[1].split('m'))
#         major = total_copy - minor
#         regions = [[chrom, start, end, major, minor]]
    return np.array(regions, dtype=[('chrom', '<U5'), ('start', int), ('end', int), ('M', int), ('m', int)])



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument('-b', '--bed', dest='bed_file', help='Path of called VCF for normal cells')
    target_group.add_argument('-r', '--region-copy', dest='region', nargs=2, help='Target region')
    parser.add_argument('-n', '--normal-vcf', dest='normal_vcf', required=True, help='Path of called VCF for normal cells')
    parser.add_argument('-t', '--tumor-vcf', dest='tumor_vcf', required=True, help='Path of called VCF for tumor cells')
    parser.add_argument('-v', '--min-num-het-variants', dest='min_num_het_variants', required=False, type=int, default=30, help='Minimum number of heterozygeous SNVs in a segment. Default: 30')
    # parser.add_argument('-c', '--min-abnormal-copy', dest='min_abnormal_copy', required=False, default='9m4', help='Minimum copy that is accepted as an abnormal segment, necessary for the allele imbalance threshold. Default: 9m4')
    parser.add_argument('-p', '--purity', dest='purity', required=False, type=float, default=1, help='Tumor purity. Default: 1')
    parser.add_argument('-S', '--sample-normal', dest='sample_normal', required=True, help='Sample normal')
    parser.add_argument('-s', '--sample-tumor', dest='sample_tumor', required=True, help='Sample tumor')
    parser.add_argument('-o', '--out', dest='out_vcf', required=True, help='Phased output VCF file path')
    parser.add_argument('-c', '--check-file', dest='out_check', required=True, help='Check results')
    parser.add_argument('-m', '--uniq-map-file', dest='uniq_map', required=True, help='Unique map regions')
    args = parser.parse_args()

#     normal_vcf = VariantFile(args.normal_vcf)
#     tumor_vcf = VariantFile(args.tumor_vcf)
#     tumor_vcf.header.formats.add('PB', '1', 'Float', 'Phred score of probability that the genotype is a balanced outcome.')
#     phased_vcf = VariantFile(args.out_vcf, 'w', header=tumor_vcf.header)
    regions = region_parser(args)
    # t, m = map(int, args.min_abnormal_copy.split('m'))

#     tumor_dna_frac = cal_tumor_dna_frac(regions, args.purity)
#     orig_ads, de_ads, n_snps = get_ad(normal_vcf, tumor_vcf, regions, tumor_dna_frac, args.sample_normal, args.sample_tumor)
#     binom_ais = cal_binom_ais_mean_all(orig_ads)
#     avg_depths = cal_avg_depths(de_ads)
#     stds = cal_ad_diff_std(orig_ads)
#     df = pd.concat(
#         [pd.DataFrame(regions, columns=['chrom', 'start', 'end', 'M', 'm']), 
#          pd.DataFrame({'n_snps' : n_snps, 
#                        'binom_ais' : binom_ais, 
#     #                    'st_binom_ais' : prep.StandardScaler()
#     #                        .fit_transform(binom_ais.reshape(-1, 1)).reshape(1, -1)[0],
#                        'avg_depths' : avg_depths, 
#     #                    'st_avg_depths' : prep.StandardScaler()
#     #                        .fit_transform(avg_depths.reshape(-1, 1)).reshape(1, -1)[0]
#                        'ad_dff_stds' : stds,
#                        'aisXstds' : binom_ais * stds,
#                       })
#         ], axis=1)
    
#     dp_ms, ais_ms, res = cnahap.check_balance(normal_vcf, tumor_vcf, args.min_num_het_variants, df)

    tumor_major_frac, tumor_minor_frac = cal_tumor_dna_frac(regions, args.purity)
#     normal_recs = np.array(list(map(get_region_recs, [normal_vcf] * regions.size, regions, [args.sample_normal] * regions.size)))
#     tumor_recs = np.array(list(map(get_region_recs, [tumor_vcf] * regions.size, regions, [args.sample_tumor] * regions.size)))
#     tumor_recs_used = np.array(list(map(filter_region_recs, normal_recs, tumor_recs)))
    normal_recs = get_recs(args.normal_vcf, args.uniq_map)
    print('Normal rec done')
    tumor_recs = get_recs(args.tumor_vcf, args.uniq_map)
    print('Tumor rec done')
    tumor_recs_used = filter_recs(normal_recs, tumor_recs)
    tumor_recs_used_reg = np.array([get_region_recs(tumor_recs_used, region) for region in regions])
    print('Region rec done')
    de_ads = np.array(list(map(cal_region_ads_by_frac, tumor_recs_used_reg, tumor_major_frac, tumor_minor_frac)))
    binom_ais = np.array([cal_region_binom_ais_mean(ads) for ads in de_ads])
    n_snps = np.array(list(map(lambda x: x.size, de_ads)))
    avg_depths = np.array(list(map(cal_region_avg_depth, de_ads)))
    ad_diff_stds = np.array(list(map(cal_region_ad_diff_std, de_ads)))
    df = pd.concat(
        [pd.DataFrame(regions, columns=['chrom', 'start', 'end', 'M', 'm']), 
         pd.DataFrame({'n_snps' : n_snps, 
                       'binom_ais' : binom_ais, 
    #                    'st_binom_ais' : prep.StandardScaler()
    #                        .fit_transform(binom_ais.reshape(-1, 1)).reshape(1, -1)[0],
                       'avg_depths' : avg_depths, 
    #                    'st_avg_depths' : prep.StandardScaler()
    #                        .fit_transform(avg_depths.reshape(-1, 1)).reshape(1, -1)[0]
                       'ad_diff_stds' : ad_diff_stds,
                       'aisXstds' : binom_ais * ad_diff_stds,
                      })
        ], axis=1)
    print('Stat done')
    dp_ms, ais_ms, aisXstds_ms, res = check_balance(df, args.min_num_het_variants)
    res = res.sort_values(['chrom', 'start'])
    tumor_recs_used_reg = tumor_recs_used_reg[res.index]
    de_ads = de_ads[res.index]
    res.to_csv(args.out_check, index=False, sep='\t')
#     will_phase_df = check_balance(normal_vcf, tumor_vcf, args.min_num_het_variants, info_df)
#     will_phase_df.to_csv(args.out_check, index=False, sep='\t')
    print('Segment check done')
    
    ba = res.loc[lambda row: row.n_snps > args.min_num_het_variants].loc[lambda row: row.m > 0].loc[lambda row: row.M == row.m]
    imba = res.loc[lambda row: row.n_snps > args.min_num_het_variants].loc[lambda row: row.m > 0].loc[lambda row: row.M != row.m]
    print(f'TP: {imba.phased.sum()} / {len(imba)} = {imba.phased.sum() / len(imba)}')
    print(f'TN: {(ba.phased == False).sum()} / {len(ba)} = {(ba.phased == False).sum() / len(ba)}')
    print(f'FP: {ba.phased.sum()} / {len(ba)} = {ba.phased.sum() / len(ba)}')
    print(f'FN: {(imba.phased == False).sum()} / {len(imba)} = {(imba.phased == False).sum() / len(imba)}')
    print(f'ACC: {imba.phased.sum() + (ba.phased == False).sum()} / {len(imba) + len(ba)} = {(imba.phased.sum() + (ba.phased == False).sum()) / (len(imba) + len(ba))}')
    print(f'PRE: {imba.phased.sum()} / {imba.phased.sum() + ba.phased.sum()} = {imba.phased.sum() / (imba.phased.sum() + ba.phased.sum())}')
    
    phase_segment(args.out_vcf, res, tumor_recs_used_reg, de_ads)
    
#     for i in range(len(regions)):
#         print(f'{regions[i][0]} {regions[i][1]} {regions[i][2]} {regions[i][3]} {regions[i][4]} {len(ads[i])} {binom_ais[i] if len(ads[i]) > 0 else "NA"} {avg_depths[i]} {will_phase[i]}')
        # print(f'{regions[i][0]} {regions[i][1]} {regions[i][2]} {len(ads[i])} {binom_ais[i] if len(ads[i]) > 0 else "NA"} {avg_depths[i]} {will_phase[i]}')
#     for i, region in enumerate(regions):
#         if not will_phase[i]:
#             continue
#         phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, args.purity, args.sample_normal, args.sample_tumor)
#     for row in will_phase_df.itertuples():
#         if row.phased:
#             region = (row.chrom, row.start, row.end)
#             phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, args.purity, args.sample_normal, args.sample_tumor)
#     phased_vcf.close()
# 