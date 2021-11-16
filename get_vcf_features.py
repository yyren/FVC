from __future__ import (absolute_import, division, print_function,
                        unicode_literals, with_statement)
import argparse
import numpy as np
import collections
import math
import re
import itertools
import datetime
import sys


class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def get_idx_value(value,type):
    values=value.strip().split(',')
    if type == 'first':
        return(values[0])
    if type == 'second':
        return(values[1])
    elif type == 'freq':
        total = float(values[0])+int(values[1])
        if total == 0:
            freq=0
        else:
            freq=('%.3f' %(float(values[1])/total))
        return(freq)
    elif type == 'str_len':
        return(len(value))

def get_AF_adj(AF,AD):
    AD=float(AD)
    AF=float(AF)
    if AD <= 1:
        return(0)
    else:
        AF_adj=('%.3f' %(float((AD-2)/(AD*(1+(math.exp(-AF)))))))
        return(AF_adj)


#MBQ=MBQ[1] MFRL[1] RPA1=RPA[0] RPA2=RPA[1] AU: length ClusterEvent 
# af_idx=len(Infor_feature)
# col_numb=af_idx+1
def get_annotations(QUAL,INFO,FORMAT,Value,ALT,col_numb,af_idx,Infor_feature):
    idx='first'
    if (re.findall(r'\*',ALT)):
        alts=ALT.strip().split(',')
        if alts[0] == '*':
            idx='second'
    annotation_feature=np.zeros((1,col_numb))
    infors=INFO.strip().split(';')
    for i in range(len(infors)):
        if re.findall(r'=',infors[i]):
            infor_recs=infors[i].strip().split('=')
            if infor_recs[0] in Infor_feature:
                if infor_recs[0] == 'RPA':
                    RPA_a1=get_idx_value(infor_recs[1],'first')
                    RPA_a2=get_idx_value(infor_recs[1],'second')
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=RPA_a1
                    annotation_feature[0][Infor_feature[infor_recs[0]]+1]=RPA_a2
                elif infor_recs[0] == 'RU':
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=get_idx_value(infor_recs[1],'str_len')
                elif infor_recs[0] == 'MBQ' or infor_recs[0] == 'MFRL':
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=get_idx_value(infor_recs[1],'second')
                elif  infor_recs[0] == 'MPOS':
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=get_idx_value(infor_recs[1],idx)
                elif infor_recs[0] == 'AS_SB_TABLE':
                    values=infor_recs[1].split('|')
                    alt_strand_covs=values[1].split(',')
                    max_value=max(float(alt_strand_covs[0]),float(alt_strand_covs[1]))
                    min_value=min(float(alt_strand_covs[0]),float(alt_strand_covs[1]))
                    strand_bias=0.0
                    if max_value != 0:
                        strand_bias=('%.3f' %(float(min_value)/float(max_value)))
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=strand_bias
                elif infor_recs[0] == 'AS_UNIQ_ALT_READ_COUNT':
                    ad_cov=infor_recs[1]
                    if re.findall(r'|',infor_recs[1]):
                        ad_covs=infor_recs[1].split('|')
                        ad_cov=ad_covs[0]
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=ad_cov
                else:
                    annotation_feature[0][Infor_feature[infor_recs[0]]]=infor_recs[1]
                
    formats=FORMAT.strip().split(':')
    values=Value.strip().split(':')
    for i in range(len(formats)):
        if formats[i] == 'AD':
            if (not (re.findall(r'\,',values[i]))):
                print('warning:'+Value)
                values[i]='0,'+str(values[i])
            AF=get_idx_value(values[i],'freq')
            AD=get_idx_value(values[i],'second')
            AF_adj=get_AF_adj(AF,AD)
            annotation_feature[0][af_idx]=AF_adj
    return annotation_feature

def out_record(vcf_file,out_tensor,tag='0',model='gatk',max_sample_numb=99999999999):
    begin2end = MagicDict()
    sample_numb=1
    feature = open(out_tensor,'w')
    time_stamp = datetime.datetime.now()
    print('Starting get features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
    
    Infor_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'RU':10, 'SOR':11, 'FS':12, 'AS_SB_TABLE':13, 'AS_UNIQ_ALT_READ_COUNT':14}
    if re.findall(r'gatk',model,flags=re.IGNORECASE):
        Infor_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'RU':10, 'SOR':11, 'FS':12, 'QD':13, 'ExcessHet':14, 'GQ_MEAN':15, 'AS_SB_TABLE':16, 'AS_UNIQ_ALT_READ_COUNT':17}
    elif re.findall(r'varscan',model,flags=re.IGNORECASE):
        Infor_feature = {'MQ':0, 'MQRankSum':1, 'MQ0':2, 'MBQ':3, 'BaseQRankSum':4, 'MFRL':5, 'MPOS':6, 'ReadPosRankSum':7, 'RPA':8, 'RPA2':9, 'RU':10, 'SOR':11, 'FS':12, 'GQ_MEAN':13, 'AS_SB_TABLE':14, 'AS_UNIQ_ALT_READ_COUNT':15}
    
    af_idx=len(Infor_feature)
    col_numb=af_idx+1
    
    with open(vcf_file) as f:
        for row in f.readlines():
            if not (re.findall(r'^#',row)):
                if sample_numb <= max_sample_numb:
                    row = row.strip().split()
                    chrom=re.sub('chr','',row[0])#chr1->1
                    output_line = []
                    output_line.append( "%d %s %d %d %s %s" % (int(tag),chrom,int(row[1]),int(row[1]),str(row[3]),str(row[4])) )
                    annotation_feature = get_annotations(row[5],row[7],row[8],row[9],row[4],col_numb,af_idx,Infor_feature)
                    for record in np.reshape(annotation_feature, 1*col_numb):
                        output_line.append("%.4f" % record)
                    out_line=" ".join(output_line)
                    feature.write(out_line+'\n')
                else:
                    break
                sample_numb +=1
    f.close()
    feature.close()
    time_stamp = datetime.datetime.now()
    print('Finish get features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')

def main():
    parser = argparse.ArgumentParser(description='filtering the false variants in NGS variant calling results' )
    parser.add_argument('--in_file', type=str, default="input.vcf", help="the vcf file need to be quality control, default:input.vcf")
    parser.add_argument('--out_file', type=str, default="./", help="the quality controled vcf file, default:./")
    parser.add_argument('--caller', type=str, default="gatk", help="variant caller to be selected: gatk,varscan,others, default: gatk")
    parser.add_argument('--tag', type=str, default="0", help="fp:0, tp:1, default: 0")
    
    args = parser.parse_args()
    vcf_file=args.in_file
    out_file=args.out_file
    tag=args.tag
    caller=args.caller
    out_record(vcf_file,out_file, tag, caller)

if __name__ == '__main__': main()
