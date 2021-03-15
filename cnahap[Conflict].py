import argparse
import gzip
import numpy as np
from scipy import stats
from pysam import VariantFile
import pandas as pd

# def get_major_minor(record):
#     dp4 = record.info['DP4']
#     dp_ref = sum(dp4[0:2])
#     dp_alt = sum(dp4[2:4])
#     return max([dp_ref, dp_alt]), min([dp_ref, dp_alt])

# def get_region_major_minor_diff(normal_vcf, tumor_vcf, region):
#     diff = []
#     for normal_record in normal_vcf.fetch(region[0], region[1] - 1, region[2]):
#         # include homozygotes
#         sample = normal_record.samples.keys()[0]
#         # filter out homozygotes
#         if sum(normal_record.samples[sample]['GT']) != 1:
#             continue
#         try:
#             tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
#             major, minor = get_major_minor(tumor_record)
#             diff.append(major - minor)
#         except:
#             pass
#     return diff

def get_major_minor(normal_vcf, tumor_vcf, region):
    major_l = []
    minor_l = []
    for normal_record in normal_vcf.fetch(region[0], region[1] - 1, region[2]):
        sample = normal_record.samples.keys()[0]
        # filter out homozygotes
        if sum(normal_record.samples[sample]['GT']) != 1:
            continue
        try:
            tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
            if sum(tumor_record.samples[sample]['GT']) != 1:
                continue
            # dp4 = tumor_record.info['DP4']
            # ref = sum(dp4[0:2])
            # alt = sum(dp4[2:4])
            # major_l.append(max([ref, alt]))
            # minor_l.append(min([ref, alt]))
            ad = tumor_record.samples[sample]['AD']
            # ref = sum(dp4[0:2])
            # alt = sum(dp4[2:4])
            major_l.append(max(ad))
            minor_l.append(min(ad))
        except:
            pass
    return np.array(major_l), np.array(minor_l)

def cal_rank_mean(ai, thres):
    if len(ai) == 0:
        return -1
    r = [int(i / thres) + 1 for i in ai]
    r = sorted([(i, r.count(i)) for i in set(r)], key=lambda x: x[1], reverse=True)
    w = {}
    cw = len(r) - 1
    for e in r:
    	w.update({e[0] : 2 ** cw})
    	cw -= 1
    return sum([i * w[int(i / thres) + 1] for i in ai]) / sum([w[int(i / thres) + 1] for i in ai])

def cal_rank_mean_by_ploidy(ai, max_abnormal_ploidy):
    if len(ai) == 0:
        return -1
    ai_range = [(max_abnormal_ploidy - m - m) / (max_abnormal_ploidy - m) for m in range(1, int(max_abnormal_ploidy / 2) + 1)]
    ai_range = [0] + ai_range + [1]
    s = 0
    sw = 0
    for e in ai:
        idx = 0
        while e > ai_range[idx]:
            idx += 1
        sw += idx
        s += e * idx
    return s / sw

def check_balance(normal_vcf, tumor_vcf, regions, min_num_het_variants, min_abnormal_ai):
    print('Extracting major and minor')
    major_minors = [get_major_minor(normal_vcf, tumor_vcf, region) for region in regions]
    print('Calculating AI')
    ais = [(e[0] - e[1]) / e[0] for e in major_minors]
    print('Calculating rank mean')
    rank_means = [cal_rank_mean(ai, min_abnormal_ai) for ai in ais]
    print('Calculating average depth')
    avg_depths = [(e[0] + e[1]).mean() if len(e[0]) > 0 else -1 for e in major_minors]

    are_balanced_before = [True if rank_means[i] <= min_abnormal_ai else False for i in range(len(rank_means))]
    normal_depths = [dp for is_bal, dp in zip(are_balanced_before, avg_depths) if is_bal]
    normal_avg_depth = sum(normal_depths) / len(normal_depths) / 2
    are_2copy = [True if round(dp / normal_avg_depth) == 2 else False for dp in avg_depths]
    have_enough_snv = [True if len(ai) >= min_num_het_variants else False for ai in ais]
    # are_balanced_after = None

    while True:
        normal_depths = [dp for is_2copy, dp in zip(are_2copy, avg_depths) if is_2copy]
        normal_avg_depth = sum(normal_depths) / len(normal_depths) / 2
        are_2copy_tmp = [True if round(dp / normal_avg_depth) == 2 else False for dp in avg_depths]
        if are_2copy_tmp == are_2copy:
            break
        else:
            are_2copy = are_2copy_tmp


    # while are_balanced_after != are_balanced_before:
    #     are_balanced_after = are_balanced_before
    #     normal_depths = [dp for is_bal, dp in zip(are_balanced_before, avg_depths) if is_bal]
    #     normal_avg_depth = sum(normal_depths) / len(normal_depths) / 2
    #     are_2copy = [True if round(dp / normal_avg_depth) == 2 else False for dp in avg_depths]
    
    for i in range(len(ais)):
        print(f'{regions[i][0]} {regions[i][1]} {regions[i][2]} {regions[i][3]} {regions[i][4]} {len(ais[i])} {rank_means[i] if len(ais[i]) > 0 else "NA"} {avg_depths[i]} {are_balanced_before[i]} {are_2copy[i]} {have_enough_snv[i]} {not are_balanced_before[i] and not are_2copy[i] and have_enough_snv[i]} {np.std(ais[i])}')

    return [True if not e[0] and not e[1] and e[2] else False for e in zip(are_balanced_before, are_2copy, have_enough_snv)], major_minors

def phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, major_minors, min_num_het_variants, min_abnormal_ai):
    # diff = get_region_major_minor_diff(normal_vcf, tumor_vcf, region)
    # print(f'Processing: {region[0]}:{region[1]}-{region[2]} with min_het {min_num_het_variants} and min_ai {min_abnormal_ai}')
    # major_l, minor_l = get_major_minor(normal_vcf, tumor_vcf, region)
    # diff = major_l - minor_l
    # ai = diff / major_l
    # if len(major_l) < min_num_het_variants or is_balance(ai, min_abnormal_ai):
    #     print(f'{region[0]} {region[1]} {region[2]} {region[3]} {region[4]} {len(major_l)} {cal_2p_mean(ai, min_abnormal_ai) if len(ai) > 0 else "NA"} skip.')
    #     # print(f'Skip: num_het = {len(major_l)}, mean_ai = {ai.mean()}')
    #     # Skip the segment whose number of heterozygotes is too small.
    #     # Or the mean allele imbalance is smaller than the lowest resolution for abnormal, i.e., the segment is treated as normal.
    #     return
    # print(f'{region[0]} {region[1]} {region[2]} {region[3]} {region[4]} {len(major_l)} {cal_2p_mean(ai, min_abnormal_ai) if len(ai) > 0 else "NA"} continue.')
    
    diff_thres = np.quantile(major_minors[0] - major_minors[1], 0.05)
    # if region[3] + region[4] != 2:
    # return
    segment_first = True
    num = 0
    for normal_record in normal_vcf.fetch(region[0], region[1] - 1, region[2]):
        sample = normal_record.samples.keys()[0]
        if sum(normal_record.samples[sample]['GT']) != 1:
            # filter out homozygotes
            continue
        if len(normal_record.ref) > 1 or any(map(lambda x : len(x) > 1, normal_record.alts)):
            # Skip the non-SNV loci
            continue
        # catch exception that the record doesn't exist
        try:
            tumor_record = next(tumor_vcf.fetch(normal_record.chrom, normal_record.pos - 1, normal_record.pos))
            # dp4 = tumor_record.info['DP4']
            # dp_ref = sum(dp4[0:2])
            # dp_alt = sum(dp4[2:4])
            # major, minor = region[3:5]
            # if dp_ref + dp_alt < 5:
            #     continue
            # if (max([dp_ref, dp_alt]) - min(dp_ref, dp_alt)) / max([dp_ref, dp_alt]) < (max([major, minor]) - min(major, minor)) / max([major, minor]):
            #     print(f'{major} {minor} {(major - minor) / major} {dp_ref} {dp_alt} {(max([dp_ref, dp_alt]) - min(dp_ref, dp_alt)) / max([dp_ref, dp_alt])}')
            #     continue
            # if max([dp_ref, dp_alt]) - min([dp_ref, dp_alt]) < diff_thres:
            #     continue
            # if dp_ref > dp_alt:
            #     # hap = (r.chrom, r.pos, 0, 1)
            #     tumor_record.samples[sample]['GT'] = (0, 1)
            # else:
            #     # hap = (r.chrom, r.pos, 1, 0)
            #     tumor_record.samples[sample]['GT'] = (1, 0)
            ad = tumor_record.samples[sample]['AD']
            if sum(ad) < 5:
                continue
            if max(ad) - min(ad) < diff_thres:
                continue
            if ad[0] > ad[1]:
                # hap = (r.chrom, r.pos, 0, 1)
                tumor_record.samples[sample]['GT'] = (0, 1)
            else:
                # hap = (r.chrom, r.pos, 1, 0)
                tumor_record.samples[sample]['GT'] = (1, 0)

            if segment_first:
                tumor_record.samples[sample].phased = False
            else:
                tumor_record.samples[sample].phased = True
            phased_vcf.write(tumor_record)
            num += 1
            segment_first = False
        except:
            pass
    print(f'{region[0]} {region[1]} {region[2]} {region[3]} {region[4]} {len(major_minors[0])} {num}')
        
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
            regions.append([chrom, start, end, major, minor])
            line = fin.readline()[:-1]
    else:
        chrom = args.region[0].split(':')[0]
        start, end = map(int, args.region.split(':')[1].split('-'))
        total_copy, minor = map(int, args.region[1].split('m'))
        major = total_copy - minor
        regions = [[chrom, start, end, major, minor]]
    return regions
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument('-b', '--bed', dest='bed_file', help='Path of called VCF for normal cells')
    target_group.add_argument('-r', '--region-copy', dest='region', nargs=2, help='Target region')
    parser.add_argument('-n', '--normal-vcf', dest='normal_vcf', required=True, help='Path of called VCF for normal cells')
    parser.add_argument('-t', '--tumor-vcf', dest='tumor_vcf', required=True, help='Path of called VCF for tumor cells')
    parser.add_argument('-v', '--min-num-het-variants', dest='min_num_het_variants', required=False, type=int, default=30, help='Minimum number of heterozygeous SNVs in a segment. Default: 30')
    parser.add_argument('-c', '--min-abnormal-copy', dest='min_abnormal_copy', required=False, default='9m4', help='Minimum copy that is accepted as an abnormal segment, necessary for the allele imbalance threshold. Default: 9m4')
    parser.add_argument('-o', '--out', dest='out_vcf', required=True, help='Phased output VCF file path')
    args = parser.parse_args()

    normal_vcf = VariantFile(args.normal_vcf)
    tumor_vcf = VariantFile(args.tumor_vcf)
    phased_vcf = VariantFile(args.out_vcf, 'w', header=tumor_vcf.header)
    regions = region_parser(args)
    t, m = map(int, args.min_abnormal_copy.split('m'))
    will_phase, major_minors = check_balance(normal_vcf, tumor_vcf, regions, args.min_num_het_variants, (t - m - m) * 1.0 / (t - m))
    for i, region in enumerate(regions):
        if not will_phase[i]:
            continue
        phase_segment(normal_vcf, tumor_vcf, phased_vcf, region, major_minors[i], args.min_num_het_variants, (t - m - m) * 1.0 / (t - m))
    phased_vcf.close()