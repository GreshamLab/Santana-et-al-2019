# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:29:45 2019

novembro - a small script for a three-way chi-squared test on downsampled OTU frequency tables,
followed by a Mann Whiteny U test. Statistically significant taxa are output into notched boxplots.

@author: Pieter Spealman ps163@nyu.edu
"""
import argparse
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

#python novembro.py -i mangrove/Otu_frequency.csv -t mangrove/corrected_taxonomy.tsv -o mangrove/

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_otu_frequency_file_name')
parser.add_argument('-t', '--taxonomy_file_name')
parser.add_argument('-o', '--outfile_path')
args = parser.parse_args()

# in file:
otu_file = open(args.input_otu_frequency_file_name)

seco_1, seco_2, seco_3 = 0, 0, 0 
intertidal_1, intertidal_2, intertidal_3 = 0, 0, 0
submerged_1, submerged_2, submerged_3 = 0, 0, 0

for line in otu_file:
    if line[0]!='#':
        #
        seco_1 += int(line.split(',')[1])
        seco_2 += int(line.split(',')[2])
        seco_3 += int(line.split(',')[3])
        #
        intertidal_1 += int(line.split(',')[4])
        intertidal_2 += int(line.split(',')[5])  
        intertidal_3 += int(line.split(',')[6])
        #
        submerged_1 += int(line.split(',')[7])
        submerged_2 += int(line.split(',')[8])       
        submerged_3 += int(line.split(',')[9])
        
otu_file.close()

global_min = min(seco_1,seco_2,seco_3,
                 intertidal_1,intertidal_2,intertidal_3,
                 submerged_1,submerged_2,submerged_3)

#Define number of observations per site for the purposes of downsampling
seco_cor_1, seco_cor_2, seco_cor_3 = global_min/seco_1, global_min/seco_2, global_min/seco_3 
intertidal_cor_1, intertidal_cor_2, intertidal_cor_3 = global_min/intertidal_1, global_min/intertidal_2, global_min/intertidal_3
submerged_cor_1, submerged_cor_2, submerged_cor_3 = global_min/submerged_1, global_min/submerged_2, global_min/submerged_3

#store otu names
otu_set = set()
#store otu data
otu_dict = {}
otu_raw_dict = {}

otu_file = open(args.input_otu_frequency_file_name)
for line in otu_file:
    if line[0]!='#':
        otu = line.split(',')[0]
        #
        seco_1 = int(line.split(',')[1])
        seco_2 = int(line.split(',')[2])
        seco_3 = int(line.split(',')[3])
        raw_seco = [seco_1, seco_2, seco_3]
        
        intertidal_1 = int(line.split(',')[4])
        intertidal_2 = int(line.split(',')[5])
        intertidal_3 = int(line.split(',')[6])
        raw_intertidal = [intertidal_1, intertidal_2, intertidal_3]
        
        submerged_1 = int(line.split(',')[7])
        submerged_2 = int(line.split(',')[8])
        submerged_3 = int(line.split(',')[9])
        raw_submerged = [submerged_1, submerged_2, submerged_3]
        otu_raw_dict[otu] = [raw_seco, raw_intertidal, raw_submerged]

        #
        seco_1 = round(int(line.split(',')[1])*seco_cor_1)
        seco_2 = round(int(line.split(',')[2])*seco_cor_2)
        seco_3 = round(int(line.split(',')[3])*seco_cor_3)
        seco = [seco_1, seco_2, seco_3]
        #
        intertidal_1 = round(int(line.split(',')[4])*intertidal_cor_1)
        intertidal_2 = round(int(line.split(',')[5])*intertidal_cor_2)
        intertidal_3 = round(int(line.split(',')[6])*intertidal_cor_3)
        intertidal = [intertidal_1, intertidal_2, intertidal_3]
        #
        submerged_1 = round(int(line.split(',')[7])*submerged_cor_1)
        submerged_2 = round(int(line.split(',')[8])*submerged_cor_2)
        submerged_3 = round(int(line.split(',')[9])*submerged_cor_3)
        submerged = [submerged_1, submerged_2, submerged_3]
        #
        otu_set.add(otu)
        
        if otu not in otu_dict:
            otu_dict[otu] = [seco, intertidal, submerged]
        else:
            print('error: duplicate otus present', otu)
            
otu_file.close()
# 

def test_difference(mean_one, mean_two, std_one, std_two):
    if abs(mean_one-mean_two) > std_one and abs(mean_one-mean_two) > std_two:
        return(True)
    else:
        return(False)

def logic_bit(sc_beyond, ss_beyond, is_beyond):
    site = 'complicated'
    
    if sc_beyond and ss_beyond and not is_beyond:
         site = 'seco'
        
    if is_beyond and sc_beyond and not ss_beyond:
        site = 'intertidal'
            
    if ss_beyond and is_beyond and not sc_beyond:
        site = 'submerged'
    
    return(site)
        
outfile = open('outstanding_otus.tab', 'w')
header = ('#site\totu\tpval\tseco_1\tseco_2\tseco_3\tseco_mean\tseco_std\tintertidal_1\tintertidal_2\tintertidal_3\tintertidal_mean\tintertidal_std\tsubmerged_1\tsubmerged_2\tsubmerged_3\tsubmerged_mean\tsubmerged_std\n')
outfile.write(header)

def round_std_test(sample_set):
    sample_mean = np.mean(sample_set)
    sample_std = np.std(sample_set)
    
    ct = 0
    for each in sample_set:
        if abs(each-sample_mean) < sample_std:
            ct+=1
            
    if ct >= 2:
        return(True)
    else:
        return(False)
    
for otu in otu_set:
    otu_array = otu_dict[otu]
    #
    seco_mean = np.mean(otu_array[0][:3])
    seco_std = np.std(otu_array[0][:3])
    intertidal_mean = np.mean(otu_array[1][:3])
    intertidal_std = np.std(otu_array[1][:3])       
    submerged_mean = np.mean(otu_array[2][:3])
    submerged_std = np.std(otu_array[2][:3])
    
    #
    raw_otu_array = otu_raw_dict[otu]
    seco_set = raw_otu_array[0]
    intertidal_set = raw_otu_array[1]
    submerged_set = raw_otu_array[2]
    #
           
    obs = np.array([[max(seco_mean,1), (global_min - seco_mean)], [max(intertidal_mean,1), (global_min - intertidal_mean)], [max(submerged_mean,1), (global_min - submerged_mean)]])

    chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True, lambda_="log-likelihood")
    
    sc_beyond, ss_beyond, is_beyond = False, False, False

    if sum(seco_set) > 20 or sum(intertidal_set) > 20:
        _w, p_sc = stats.mannwhitneyu(seco_set, intertidal_set)
        if p_sc < 0.05:
            sc_beyond = True
        else:
            sc_beyond = False
            
    if sum(seco_set) > 20 or sum(submerged_set) > 20:
        _w, p_ss = stats.mannwhitneyu(seco_set, submerged_set)
        if p_ss < 0.05:
            ss_beyond = True
        else:
            ss_beyond = False
            
    if sum(intertidal_set) > 20 or sum(submerged_set) > 20:
        _w, p_is = stats.mannwhitneyu(intertidal_set, submerged_set)
        if p_is < 0.05:
            is_beyond = True
        else:
            is_beyond = False
            
    if (sc_beyond or ss_beyond or is_beyond) and pval <= 0.05:
        outlier_site = logic_bit(sc_beyond, ss_beyond, is_beyond)
        
        outline = ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(outlier_site, otu, pval, otu_array[0][0], otu_array[0][1], otu_array[0][2], seco_mean, seco_std, otu_array[1][0], otu_array[1][1], otu_array[1][2], intertidal_mean, intertidal_std, otu_array[2][0], otu_array[2][1], otu_array[2][2], submerged_mean, submerged_std)
        outfile.write(outline)
    
outfile.close()            
        

### Combine OTUs of the same species

taxa_file = open(args.taxonomy_file_name)

taxa_to_otu_dict = {}
otu_to_taxa_dict = {}

convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}

taxa_cutoff_name = 'species'
taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]

for line in taxa_file:
    if line[0]!='#':
        line = line.strip()
        otu = line.split('\t')[0]
        taxa = line.split('\t')[1].replace('leftNAright;','').replace('leftNAright','').replace(';','_').replace(' ','').replace('[','').replace(']','')
        if taxa[-1] == '_':
            taxa = taxa[:-1]
            
        if taxa.count('_') > taxa_cutoff_num:
            taxa_list = taxa.split('_')[:taxa_cutoff_num+1]
            taxa = ''
            for each in taxa_list:
                taxa+=str(each)+'_'
            
            if taxa[-1] == '_':
                taxa = taxa[:-1]
        
        if taxa not in taxa_to_otu_dict:
            taxa_to_otu_dict[taxa] = [otu]
        else:
            taxa_to_otu_dict[taxa].append(otu)
            
        if otu not in otu_to_taxa_dict:
            otu_to_taxa_dict[otu] = taxa
        else:
            print('err')
            
taxa_file.close()

otu_file = open(args.input_otu_frequency_file_name)

seco_1, seco_2, seco_3 = 0, 0, 0 
intertidal_1, intertidal_2, intertidal_3 = 0, 0, 0
submerged_1, submerged_2, submerged_3 = 0, 0, 0

for line in otu_file:
    if line[0]!='#':
        #
        seco_1 += int(line.split(',')[1])
        seco_2 += int(line.split(',')[2])
        seco_3 += int(line.split(',')[3])
        #
        intertidal_1 += int(line.split(',')[4])
        intertidal_2 += int(line.split(',')[5])  
        intertidal_3 += int(line.split(',')[6])
        #
        submerged_1 += int(line.split(',')[7])
        submerged_2 += int(line.split(',')[8])       
        submerged_3 += int(line.split(',')[9])
        
otu_file.close()

global_min = min(seco_1,seco_2,seco_3,
                 intertidal_1,intertidal_2,intertidal_3,
                 submerged_1,submerged_2,submerged_3)

#Define number of observations per site for the purposes of downsampling
seco_cor_1, seco_cor_2, seco_cor_3 = global_min/seco_1, global_min/seco_2, global_min/seco_3 
intertidal_cor_1, intertidal_cor_2, intertidal_cor_3 = global_min/intertidal_1, global_min/intertidal_2, global_min/intertidal_3
submerged_cor_1, submerged_cor_2, submerged_cor_3 = global_min/submerged_1, global_min/submerged_2, global_min/submerged_3

def obs_counter(list_obs_ct):
    obs = 0
    for obs_ct in list_obs_ct:
        if obs_ct > 0:
            obs+=1
    
    if obs > 1:
        return(True)
    else:
        return(False)
    
#store otu names
taxa_set = set()
#store otu data
taxa_dict = {}
taxa_dict['total']=[0, 0, 0]

taxa_raw_dict = {}

otu_file = open(args.input_otu_frequency_file_name)

for line in otu_file:
    if line[0]!='#':
        otu = line.split(',')[0]
        taxa = otu_to_taxa_dict[otu]
        #
        seco_1 = int(line.split(',')[1])
        seco_2 = int(line.split(',')[2])
        seco_3 = int(line.split(',')[3])
        raw_seco = [seco_1, seco_2, seco_3]
        
        if seco_1 > 0 and seco_1 < 1:
            print(line)
            1/0
            
        if seco_2 > 0 and seco_2 < 1:
            print(line)
            1/0
            
        if seco_3 > 0 and seco_3 < 1:
            print(line)
            1/0
         
        intertidal_1 = int(line.split(',')[4])
        intertidal_2 = int(line.split(',')[5])
        intertidal_3 = int(line.split(',')[6])
        raw_intertidal = [intertidal_1, intertidal_2, intertidal_3]
        
        submerged_1 = int(line.split(',')[7])
        submerged_2 = int(line.split(',')[8])
        submerged_3 = int(line.split(',')[9])
        raw_submerged = [submerged_1, submerged_2, submerged_3]
                
        if taxa not in taxa_raw_dict:
            taxa_raw_dict[taxa] = [raw_seco, raw_intertidal, raw_submerged]
        else:
            taxa_raw_dict[taxa][0] += raw_seco
            taxa_raw_dict[taxa][1] += raw_intertidal
            taxa_raw_dict[taxa][2] += raw_submerged

        #
        seco_1 = (int(line.split(',')[1])*seco_cor_1)
        seco_2 = (int(line.split(',')[2])*seco_cor_2)
        seco_3 = (int(line.split(',')[3])*seco_cor_3)
        seco = [seco_1, seco_2, seco_3]
        #
        intertidal_1 = (int(line.split(',')[4])*intertidal_cor_1)
        intertidal_2 = (int(line.split(',')[5])*intertidal_cor_2)
        intertidal_3 = (int(line.split(',')[6])*intertidal_cor_3)
        intertidal = [intertidal_1, intertidal_2, intertidal_3]
        #
        submerged_1 = (int(line.split(',')[7])*submerged_cor_1)
        submerged_2 = (int(line.split(',')[8])*submerged_cor_2)
        submerged_3 = (int(line.split(',')[9])*submerged_cor_3)
        submerged = [submerged_1, submerged_2, submerged_3]
        #
        taxa_set.add(taxa)
                
        if taxa not in taxa_dict:
            taxa_dict[taxa] = [seco, intertidal, submerged]
        else:
            taxa_dict[taxa][0] += seco
            taxa_dict[taxa][1] += intertidal
            taxa_dict[taxa][2] += submerged
            
        taxa_dict['total'][0] += sum(seco)
        taxa_dict['total'][1] += sum(intertidal)
        taxa_dict['total'][2] += sum(submerged)

otu_file.close()

#        
outfile_name = ('{}site_specific_{}_enrichment.tab').format(args.outfile_path, taxa_cutoff_name)
outfile = open(outfile_name, 'w')
header = ('#taxa\tuid\tchi_pval\tMWU_pval_seco_intertidal\tMWU_pval_intertidal_submerged\tMWU_pval_submerged_seco\n')
outfile.write(header)

uid = 0

max_obs = 0
figure_dict = {}

def criteria(seco_set, intertidal_set, submerged_set, p_sc, p_is, p_ss, pval, pct_effect_size=0.05):
    pass_set = []
    log_fold_diff = []
    
    if pval <= 0.05:
        effect_size = abs(np.log10((pct_effect_size*1000+1000)/1000))
        
        log_seco_median = np.log10(max(np.mean(seco_set),1))
        log_intertidal_median = np.log10(max(np.mean(intertidal_set),1))
        log_submerged_median = np.log10(max(np.mean(submerged_set),1))
            
        w_sc = (log_seco_median-log_intertidal_median)
        w_is = (log_intertidal_median-log_submerged_median)
        w_ss = (log_submerged_median-log_seco_median)
                
        log_fold_diff = [w_sc, w_is, w_ss]
        
        if (p_sc <= 0.05) and abs(w_sc) >= effect_size:
            pass_set.append('sc')
                
        if (p_is <= 0.05) and abs(w_is) >= effect_size:
            pass_set.append('is')
        
        if (p_ss <= 0.05) and abs(w_ss) >= effect_size:
            pass_set.append('ss')
            
    if len(pass_set) > 1:
        return(True, log_fold_diff, pass_set)
    else:
        return(False, log_fold_diff, pass_set)

def test_max(set_list, max_obs):
    for each_set in set_list:
        if max(each_set) >= max_obs:
            max_obs = max(each_set)
    return(max_obs)

def return_log10(each_set):
    new_set = []
    for each_obs in each_set:
        if each_obs <= 0:
            each_obs = 0
        else:
            each_obs = np.log10(each_obs)
        
        new_set.append(each_obs)
    
    return(new_set)
    
figure_dict = {}

total_seco = taxa_dict['total'][0]
total_intertidal = taxa_dict['total'][1]
total_submerged = taxa_dict['total'][2]

for taxa in taxa_set:
    
    seco_array, intertidal_array, submerged_array = taxa_dict[taxa]
    
    seco_mean = np.mean(seco_array)
    seco_std = np.std(seco_array)
    intertidal_mean = np.mean(intertidal_array)
    intertidal_std = np.std(intertidal_array)       
    submerged_mean = np.mean(submerged_array)
    submerged_std = np.std(submerged_array)
    
    raw_taxa_array = taxa_raw_dict[taxa]
    seco_set = raw_taxa_array[0]
    intertidal_set = raw_taxa_array[1]
    submerged_set = raw_taxa_array[2]
    
    seco_obs = obs_counter(seco_set)
    intertidal_obs = obs_counter(intertidal_set)
    submerged_obs = obs_counter(submerged_set)
        
    if seco_obs or intertidal_obs or submerged_obs:                
        obs = np.array([[max(seco_mean,1), (total_seco - sum(seco_array))], [max(intertidal_mean,1), (total_intertidal - sum(intertidal_array))], [max(submerged_mean,1), (total_submerged - sum(submerged_array))]])
    
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True, lambda_="log-likelihood")
                        
        if sum(seco_set) > 20 or sum(intertidal_set) > 20:
            _w, p_sc = stats.mannwhitneyu(seco_set, intertidal_set)
                
        if sum(seco_set) > 20 or sum(submerged_set) > 20:
            _w, p_ss = stats.mannwhitneyu(seco_set, submerged_set)

        if sum(intertidal_set) > 20 or sum(submerged_set) > 20:
            _w, p_is = stats.mannwhitneyu(intertidal_set, submerged_set)
            
        pass_criteria, log_fold_difference, pass_set = criteria(seco_array, intertidal_array, submerged_array, p_sc, p_is, p_ss, pval)

        if pass_criteria:
           
            log_seco_set = return_log10(seco_set)
            log_intertidal_set = return_log10(intertidal_set)
            log_submerged_set = return_log10(submerged_set)
            
            max_obs = test_max([log_seco_set,log_intertidal_set,log_submerged_set], max_obs)
            
            if taxa not in figure_dict:
                figure_dict[taxa] = [log_seco_set, log_intertidal_set, log_submerged_set]
            else:
                print('error')
                1/0

            outline = ('{}\t{}\t{}\t{}\t{}\t{}\n').format(taxa, uid, pval, p_sc, p_is, p_ss)
            outfile.write(outline)
            
            uid += 1
            
outfile.close()
            
for taxa, set_list in figure_dict.items():
    outfile_name = ('{}site_specific_{}_enrichment_{}.pdf').format(outfile_path, taxa_cutoff_name, taxa)
    
    seco_set = set_list[0]
    intertidal_set = set_list[1]
    submerged_set = set_list[2]   
    
    seco_no_zeros = [x for x in seco_set if x > 0] 
    intertidal_no_zeros = [x for x in intertidal_set if x > 0]
    submerged_no_zeros = [x for x in submerged_set if x > 0]
    
    axes = plt.gca()
    axes.set_ylim([-1, max_obs])
    
    print(taxa)
    print(len(seco_no_zeros), len(intertidal_no_zeros), len(submerged_no_zeros))
    
    seco_alpha = 1 
    if len(seco_no_zeros) > 0:
        seco_mean = np.mean(seco_no_zeros)
        seco_alpha = max((1/len(seco_no_zeros)),0.1)
        plt.scatter(1, seco_mean, c='red', marker='o')
    else:
        plt.scatter(1, 0, c='red', marker='X')
        
    plt.scatter([1]*len(seco_no_zeros),seco_no_zeros, alpha=seco_alpha)
        
    intertidal_alpha = 1
    if len(intertidal_no_zeros) > 0:
        intertidal_mean = np.mean(intertidal_no_zeros)
        intertidal_alpha = max((1/len(intertidal_no_zeros)),0.1)
        plt.scatter(2, intertidal_mean, c='red', marker='o')
    else:
        plt.scatter(2, 0, c='red', marker='X')
        
    plt.scatter([2]*len(intertidal_no_zeros),intertidal_no_zeros, alpha=intertidal_alpha)
        
    submerged_alpha = 1    
    if len(submerged_no_zeros) > 0:       
        submerged_mean = np.mean(submerged_no_zeros)
        submerged_alpha = max((1/len(submerged_no_zeros)),0.1)
        plt.scatter(3, submerged_mean, c='red', marker='o')
    else:
        plt.scatter(3, 0, c='red', marker='X')  
        
    plt.scatter([3]*len(submerged_no_zeros),submerged_no_zeros, alpha=submerged_alpha)    
        
    plt.boxplot([seco_no_zeros, intertidal_no_zeros, submerged_no_zeros],1, showfliers=False)
