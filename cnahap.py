import sys
import argparse
import gzip
import numpy as np
import scipy as sp
import pysam
import pandas as pd
import subprocess
import sklearn.cluster as cl
from Bio import bgzf
import gekko

# def get_region_recs(vcf, region, sample):
# #     print(region)
#     recs = np.array([
#         (r.pos, r.ref, r.alts[0], *r.samples[sample]['GT'], *r.samples[sample]['AD'])
#         for r in vcf.fetch(region[0], int(region[1]) - 1, int(region[2]))
#         if len(r.ref) == 1 and len(r.alts) == 1 and len(r.alts[0]) == 1
#     ])
# #     print(recs)
#     return recs

CHROM_END = {
    'chr1': 249250621,
    'chr2': 243199373,
    'chr3': 198022430,
    'chr4': 191154276,
    'chr5': 180915260,
    'chr6': 171115067,
    'chr7': 159138663,
    'chr8': 146364022,
    'chr9': 141213431,
    'chr10': 135534747,
    'chr11': 135006516,
    'chr12': 133851895,
    'chr13': 115169878,
    'chr14': 107349540,
    'chr15': 102531392,
    'chr16': 90354753,
    'chr17': 81195210,
    'chr18': 78077248,
    'chr19': 59128983,
    'chr20': 63025520,
    'chr21': 48129895,
    'chr22': 51304566,
    'chrM': 155270560
}

def segmentation(regs):
    res = []
    for c in [f'chr{i}' for i in range(1, 23)]:
        c_regs = regs[regs['chrom'] == c]
        bps = np.sort(np.hstack([c_regs['start'], c_regs['end']]))
        bps = bps[~np.isin(bps, [1, CHROM_END[c]])]
        start = 1
        for p in bps:
            res.append((c, start, p - 1))
            start = p
        res.append((c, start, CHROM_END[c]))
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
        recs.append((r[0], int(r[1]), r[3], r[4], *ad, False))
    dtypes = [
        ('chrom', '<U5'), ('pos', int), 
        ('ref', '<U1'), ('alt', '<U1'), 
        ('ad1', int), ('ad2', int),
        ('to_be_phased', bool)
    ]
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
    
def mark_recs(normal_recs, tumor_recs):
    
    # Intersection between normal and tumor
    fields = ['chrom', 'pos', 'ref', 'alt']
    inter, normal_inter_idx, tumor_inter_idx = np.intersect1d(normal_recs[fields], tumor_recs[fields], return_indices=True)
    normal_recs = normal_recs[normal_inter_idx]
#     tumor_recs = tumor_recs[tumor_inter_idx]
    
    tmp = np.full(len(tumor_recs), False)
    tmp[tumor_inter_idx] = True
    
    # Mark sites to be phased
#     tumor_recs[np.logical_and(tumor_recs['ad1'] > 0, tumor_recs['ad2'] > 0)]['to_be_phased'] = True
    tmp[np.logical_or(tumor_recs['ad1'] == 0, tumor_recs['ad2'] == 0)] = False
    tumor_recs['to_be_phased'] = tmp
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

## may be some bugsd
def cal_region_ads_by_frac(recs, major_frac, minor_frac):
    if recs.size == 0 or np.isnan(major_frac):
        return recs
    ad1_ge_ad2 = recs['ad1'] >= recs['ad2']
    tmp_ad1 = np.zeros(len(recs))
    tmp_ad2 = np.zeros(len(recs))
    tmp_ad1[ad1_ge_ad2] = (recs[ad1_ge_ad2]['ad1'] * major_frac).round()
    tmp_ad2[ad1_ge_ad2] = (recs[ad1_ge_ad2]['ad2'] * minor_frac).round()
    tmp_ad1[~ad1_ge_ad2] = (recs[~ad1_ge_ad2]['ad1'] * minor_frac).round()
    tmp_ad2[~ad1_ge_ad2] = (recs[~ad1_ge_ad2]['ad2'] * major_frac).round()
    
    recs['ad1'] = tmp_ad1
    recs['ad2'] = tmp_ad2
    return recs
#     ads = get_region_ads(recs)
#     ads_major = ads.max(axis=1)
#     ads_minor = ads.min(axis=1)
#     res_ads = np.array(list(zip((ads_major * major_frac).round(), (ads_minor * minor_frac).round())))
#     return res_ads[(res_ads > 0).sum(axis=1) == 2]

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
                
def phase_segment_impute2_fmt(out_dir, samplename, res, tumor_recs_used_reg):
    with open(f'{out_dir}/{samplename}.samples', 'w') as samp_out:
        samp_out.write('sample population group sex\n')
        samp_out.write(f'{samplename} {samplename} {samplename} 2\n')
    
    total_snps = 0
    total_regs = 0
    for chrom in [f'chr{i}' for i in range(1, 23)]:
        index_chrom = res[res.phased].loc[lambda r: r.chrom == chrom].index
        recs_reg_for_phase = tumor_recs_used_reg[index_chrom]
             
        n_regs = len(recs_reg_for_phase)
        c_reg = 1
        with gzip.open(f'{out_dir}/{samplename}.{chrom}.hap.gz', 'wb') as hap_out, gzip.open(f'{out_dir}/{samplename}.{chrom}.legend.gz', 'wb') as leg_out:
            leg_out.write(b'id position a0 a1\n')
            for recs in recs_reg_for_phase:
                segment_first = True
                first_max = -1
                n_snps = len(recs)
                c_snp = 1
                for rec in recs:
                    total_snps += 1
                    print(f'Writing {c_snp}/{n_snps} SNVs in {c_reg}/{n_regs} regions in {chrom}. Total SNPs written: {total_snps}, total regions written: {total_regs} {" " * 20}\r', end='')
                    leg_out.write(f'{rec[0]}:{rec[1]}_{rec[2]}_{rec[3]} {rec[1]} {rec[2]} {rec[3]}\n'.encode())
                    if first_max < 0:
                        if rec[4] > rec[5]:
                            first_max = 0
                        else:
                            first_max = 1
                    if segment_first:
                        if first_max == 0:
                            if rec[4] > rec[5]:
                                line = '0* 1*'
                            else:
                                line = '1* 0*'
                        else:
                            if rec[4] > rec[5]:
                                line = '1* 0*'
                            else:
                                line = '0* 1*'
                        segment_first = False
                    else:
                        if first_max == 0:
                            if rec[4] > rec[5]:
                                line = '0 1'
                            else:
                                line = '1 0'
                        else:
                            if rec[4] > rec[5]:
                                line = '1 0'
                            else:
                                line = '0 1'
                    hap_out.write(f'{line}\n'.encode())
                    c_snp += 1
                c_reg += 1
                total_regs += 1

                
def phase_segment_vcf_fmt(out_fn, samplename, res, recs_used_reg):
    out_vcf = pysam.VariantFile(out_fn, 'w')
    
    # Create FILTER
    
    # Create Formats
    formats = out_vcf.header.formats
    formats.add(id='GT', number=1, type='String', description='Genotype \b')
    formats.add(id='AD', number='.', type='Integer', description='Allelic depths for the ref and alt alleles in the order listed')
    formats.add(id='DP', number=1, type='Integer', description='Approximate read depth')
    formats.add(id='PS', number=1, type='Integer', description='Phase set indicator for haplotype evaluation')
    
    # Create INFO
    info = out_vcf.header.info
    info.add(id='AC', number='A', type='Integer', description='Allele count in genotypes, for each ALT allele, in the same order as listed')
    info.add(id='AF', number='A', type='Float', description='Allele Frequency, for each ALT allele, in the same order as listed')
    info.add(id='AN', number=1, type='Integer', description='Total number of alleles in called genotypes')
    
    # Create contigs
    contigs = out_vcf.header.contigs
    for contig, length in CHROM_END.items():
        contigs.add(id=contig, length=length, assembly='hg19')
        
    # Create sample
    out_vcf.header.add_sample(samplename)
    
    total_snps = sum(map(len, recs_used_reg))
    total_regs = len(res)
    total_snps_written = 0
    ps = 1
    c_reg = 1
    for reg in res.itertuples():
        c_rec = 1
        segment_first = True
        first_max = -1
        n_snps = len(recs_used_reg[reg.Index])
        for rec in recs_used_reg[reg.Index]:
            r = out_vcf.new_record()
            r.contig = rec[0]
            r.pos = rec[1]
            r.ref = rec[2]
            r.alts = list(rec[3])
            r.info.update({
                'AC': 2,
                'AF': 0.5,
                'AN': 2
            })
            r_sample = r.samples.get(samplename)
            if not reg.phased:
                r_sample.update({
                    'GT' : (0, 1),
                    'AD' : (int(rec[4]), int(rec[5])),
                    'DP' : int(rec[4] + rec[5]),
                })
                r_sample.phased = False
            else:
                # Not phase 0 AD SNVs and non-inter SNVs
                if rec[4] == 0 or rec[5] == 0 or rec[6] == False:
                    r_sample.update({
                        'GT' : (0, 1),
                        'AD' : (int(rec[4]), int(rec[5])),
                        'DP' : int(rec[4] + rec[5]),
                    })
                    r_sample.phased = False
                else:
                    if segment_first:
                        segment_first = False
                        if rec[4] > rec[5]:
                            first_max = 0
                        else:
                            first_max = 1
                        r_sample.update({
                            'GT' : (0, 1),
                            'AD' : (int(rec[4]), int(rec[5])),
                            'DP' : int(rec[4] + rec[5]),
                            'PS' : ps
                        })
                        r_sample.phased = True
                    else:
                        if first_max == 0:
                            if rec[4] > rec[5]:
                                r_sample.update({
                                    'GT' : (0, 1),
                                    'AD' : (int(rec[4]), int(rec[5])),
                                    'DP' : int(rec[4] + rec[5]),
                                    'PS' : ps
                                })
                            else:
                                r_sample.update({
                                    'GT' : (1, 0),
                                    'AD' : (int(rec[4]), int(rec[5])),
                                    'DP' : int(rec[4] + rec[5]),
                                    'PS' : ps
                                })
                        else:
                            if rec[4] > rec[5]:
                                r_sample.update({
                                    'GT' : (1, 0),
                                    'AD' : (int(rec[4]), int(rec[5])),
                                    'DP' : int(rec[4] + rec[5]),
                                    'PS' : ps
                                })
                            else:
                                r_sample.update({
                                    'GT' : (0, 1),
                                    'AD' : (int(rec[4]), int(rec[5])),
                                    'DP' : int(rec[4] + rec[5]),
                                    'PS' : ps
                                })
                        r_sample.phased = True
            out_vcf.write(r)
            c_rec += 1
            total_snps_written += 1
            print(f'Writing {c_rec}/{n_snps} SNVs in {c_reg}/{total_regs} regions. Total SNPs written: {total_snps_written}/{total_snps}.{" " * 20}\r', end='')
        c_reg += 1
        if not segment_first:
            ps += 1
    out_vcf.close()
        
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
    parser.add_argument('-f', '--cnv-file', dest='cnv_file', required=True, help='CNV segment file from Accurity')
    parser.add_argument('-n', '--normal-vcf', dest='normal_vcf', required=True, help='Path of called VCF for normal cells')
    parser.add_argument('-t', '--tumor-vcf', dest='tumor_vcf', required=True, help='Path of called VCF for tumor cells')
    parser.add_argument('-v', '--min-num-het-variants', dest='min_num_het_variants', required=False, type=int, default=30, help='Minimum number of heterozygeous SNVs in a segment. Default: 30')
    # parser.add_argument('-c', '--min-abnormal-copy', dest='min_abnormal_copy', required=False, default='9m4', help='Minimum copy that is accepted as an abnormal segment, necessary for the allele imbalance threshold. Default: 9m4')
    parser.add_argument('-p', '--purity', dest='purity', required=True, type=float, help='Tumor purity. Default: 1')
#     parser.add_argument('-S', '--sample-normal', dest='sample_normal', required=True, help='Sample normal')
#     parser.add_argument('-s', '--sample-tumor', dest='sample_tumor', required=True, help='Sample tumor')
    parser.add_argument('-o', '--out-vcf', dest='out_fn', required=True, help='Phased output file path')
    parser.add_argument('-s', '--sample', dest='samplename', required=True, help='Sample name')
    parser.add_argument('-c', '--check-file', dest='out_check', required=True, help='Check results/Copy number estimate results')
    args = parser.parse_args()

    cnv = pd.read_csv(args.cnv_file, sep='\t')
    cnv = cnv.loc[lambda row : row.cp - np.floor(row.cp) == 0]
    cnv.chr = 'chr' + cnv.chr

    regs = []
    for row in cnv.itertuples():
        regs.append((row.chr, int(row.start), int(row.end)))

    regs = np.array(regs, dtype=[('chrom', '<U5'), ('start', int), ('end', int)])
    regs = segmentation(regs)
    
    accurity_cp = np.full(len(regs), np.nan)
    accurity_major_cp = np.full(len(regs), np.nan)
    for i, reg in enumerate(regs):
        found = cnv.loc[lambda r: r.chr == reg[0]].loc[lambda r: r.start == reg[1]]
        if len(found) > 0:
            accurity_cp[i] = found.cp.values[0]
            accurity_major_cp[i] = found.major_allele_cp.values[0]

    print('Read records')
    normal_recs = get_recs(args.normal_vcf, regs)
    tumor_recs = get_recs(args.tumor_vcf, regs)
    recs_used = filter_recs(normal_recs, tumor_recs)
    recs_used = np.array([get_region_recs(recs_used, region) for region in regs])

    print('Calculate statistics')
#     de_ads = np.array(list(map(cal_region_ads_by_frac, recs_used, np.ones(len(regs)), np.ones(len(regs)))))
    de_ads = np.array(list(map(lambda arr: np.array(arr[['ad1', 'ad2']].tolist()), recs_used)))
    binom_ais = np.array([cal_region_binom_ais_mean(ads) for ads in de_ads])
    n_snps = np.array(list(map(lambda x: x.shape[0], de_ads)))
    avg_depths = np.array(list(map(cal_region_avg_depth, de_ads)))
    M_depths = np.array(list(map(cal_region_avg_M, de_ads)))
    m_depths = np.array(list(map(cal_region_avg_m, de_ads)))
    ad_diff = np.array(list(map(cal_region_ad_diff, de_ads)))
  
    df = pd.concat(
        [pd.DataFrame(regs, columns=['chrom', 'start', 'end']), 
         pd.DataFrame({
             'accurity_cp' : accurity_cp,
             'accurity_major_cp' : accurity_major_cp,
             'n_snps' : n_snps, 
             'binom_ais' : binom_ais, 
             'avg_depths' : avg_depths, 
             'M_depths': M_depths,
             'm_depths': m_depths,
             'ad_diff' : ad_diff,
         })
        ], axis=1)
    
    print('Check normal')
    dp_ms, res = check_normal(df, args.min_num_het_variants)
    df_lt30 = res[(res.n_snps > args.min_num_het_variants)]
    hap_dp = df_lt30[(~df_lt30.dp_phased)].avg_depths.mean() / 2
    df_lt30_abnormal = df_lt30[df_lt30.dp_phased]

    print('ILP')
    # hap_dp = (df_normal.avg_depths * df_normal.n_snps).sum() / df_normal.n_snps.sum() / 2
    mod = gekko.GEKKO(remote=False) # Initialize gekko
    mod.options.SOLVER = 1  # APOPT is an MINLP solver
    mod.solver_options = ['minlp_maximum_iterations 500', \
                            # minlp iterations with integer solution
                            'minlp_max_iter_with_int_sol 200', \
                            # treat minlp as nlp
                            'minlp_as_nlp 0', \
                            # nlp sub-problem max iterations
                            'nlp_maximum_iterations 50', \
                            # 1 = depth first, 2 = breadth first
                            'minlp_branch_method 1', \
                            # maximum deviation from whole number
                            'minlp_integer_tol 0.05', \
                            # covergence tolerance
                            'minlp_gap_tol 0.01']
    M = np.full(len(res), mod.Const(0))
    # M[df_lt30_abnormal.index] = mod.Var(lb=1, integer=True)
    m = np.full(len(res), mod.Const(0))
    # m[df_lt30_abnormal.index] = mod.Var(lb=1, integer=True)
    # eM = np.full(len(res), mod.Const(0))
    # eM[df_lt30_abnormal.index] = mod.Var(lb=0)
    # em = np.full(len(res), mod.Const(0))
    # em[df_lt30_abnormal.index] = mod.Var(lb=0)

    #     p = mod.Var(value=0.1, lb=0, ub=1)
    p = args.purity
    s = []
    for row in df_lt30_abnormal.itertuples():
        M_ = mod.Var(lb=1, integer=True)
        m_ = mod.Var(lb=0, integer=True)
        eM_ = mod.Var(lb=0)
        em_ = mod.Var(lb=0)
        eD_ = mod.Var(lb=0)
        ed_ = mod.Var(lb=0)
        eB_ = mod.Var(lb=0)
    #     ea_ = mod.Var(lb=0)
        M[row.Index] = M_
        m[row.Index] = m_
    #     eM[row.Index] = eM_
    #     em[row.Index] = em_
        mod.Equation(p * hap_dp * M_ + (1 - p) * hap_dp - row.M_depths >= -eM_)
        mod.Equation(p * hap_dp * M_ + (1 - p) * hap_dp - row.M_depths <= eM_)
        mod.Equation(p * hap_dp * m_ + (1 - p) * hap_dp - row.m_depths >= -em_)
        mod.Equation(p * hap_dp * m_ + (1 - p) * hap_dp - row.m_depths <= em_)
        mod.Equation(p * (M_ + m_) * hap_dp + (1 - p) * 2 * hap_dp - row.avg_depths >= -eD_)
        mod.Equation(p * (M_ + m_) * hap_dp + (1 - p) * 2 * hap_dp - row.avg_depths  <= eD_)
        mod.Equation(p * M_ * hap_dp - p * m_ * hap_dp - row.ad_diff >= -ed_)
        mod.Equation(p * M_ * hap_dp - p * m_ * hap_dp - row.ad_diff <= ed_)
        mod.Equation((p * M_ - p * m_) / (p * M_ + 1 - p) - row.binom_ais >= -eB_)
        mod.Equation((p * M_ - p * m_) / (p * M_ + 1 - p) - row.binom_ais <= eB_)
        mod.Equation(M_ >= m_)

        s.append(np.log2(row.n_snps) * (eM_ + em_ + eD_ + ed_ + eB_))
    mod.Obj(mod.sum(s))

    mod.solve(disp=False)

    for i in range(len(M)):
        if isinstance(M[i], gekko.gk_variable.GKVariable):
            M[i] = M[i].value[0]
        else:
            M[i] = np.nan
        if isinstance(m[i], gekko.gk_variable.GKVariable):
            m[i] = m[i].value[0]
        else:
            m[i] = np.nan

    res_all = pd.concat([
        res,
        pd.DataFrame({
            'cp': M + m,
            'M': M,
            'm': m,
            'phased': np.logical_and(M > m, m > 0),
        })
    ], axis=1)

    tumor_major_frac, tumor_minor_frac = cal_tumor_dna_frac(res_all, p)
#     de_ads_est = np.array(list(map(cal_region_ads_by_frac, recs_used, tumor_major_frac, tumor_minor_frac)))

    res_all.to_csv(args.out_check, sep='\t', header=True, index=False, na_rep='na')
    
    print('Phasing')
    phase_segment_impute2_fmt(args.out_fn, args.samplename, res_all, recs_used)
