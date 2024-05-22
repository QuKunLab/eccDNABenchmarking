import os
import sys
import pandas as pd
import numpy as np
import subprocess as sp
import multiprocessing as mp
from importlib.resources import files
from os.path import exists,join,split
from .utils import utilities

class statistics(object):
    '''
    simple_ratio: [0-1]
    simple_template: bed file of simple eccDNA in bed format 'chr start end length id'
    chimeric_template: bed file of chimeric eccDNA in bed format 'chr start end length id'
    ont_mean: average of ONT reads length
    ont_sd: standard deviation of ONT reads length 
    sr_mean: average of NGS reads length
    sr_sd: standard deviation of NGS reads length
    '''
    
    def __init__(self,
                 simple_ratio = None,
                 simple_template = None,
                 chimeric_template = None,
                 ont_mean = 3000,
                 ont_std = 2500,
                 sr_mean = 400,
                 sr_std = 125,
                ):
        
        self.simple_ratio = simple_ratio
        self.simple_template = simple_template
        self.chimeric_template = chimeric_template
        self.ont_mean = ont_mean
        self.ont_std = ont_std
        self.sr_mean = sr_mean
        self.sr_std = sr_std
               
        if not self.simple_ratio:    
            if not self.simple_template and not self.chimeric_template:
                self.simple_ratio = 0.05
                self.simple_weight = self.cal_simple_weight()
                self.chimeric_weight = self.cal_chimeric_weight()
                
            elif not self.simple_template:
                self.simple_ratio = 0
                self.chimeric_weight = self.cal_chimeric_weight()
                
            elif not self.chimeric_template:
                self.simple_ratio = 1
                self.simple_weight = self.cal_simple_weight()
                
            else: 
                self.simple_ratio=self.cal_simple_ratio()
                self.simple_weight = self.cal_simple_weight()
                self.chimeric_weight = self.cal_chimeric_weight()
        else:
            self.simple_weight = self.cal_simple_weight()
            self.chimeric_weight = self.cal_chimeric_weight()
        
    def cal_simple_ratio(self):
        '''
        Function to calculate the ratio of simpleDNA molecules from template
        '''
        simple = pd.read_csv(self.simple_template,header=None,sep='\t')
        chimeric = pd.read_csv(self.chimeric_template,header=None,sep='\t')
        simple_count = len(simple)
        chimeric_count = len(chimeric[4].unique())
        return simple_count/(simple_count + chimeric_count)
    
    def cal_simple_weight(self):
        '''
        Function to calculate the chr and length weight of simple DNA molecules from template (if not ).
        '''
        if not self.simple_template:
            self.simple_template = join(split(__file__)[0], 'resource', 'template', 'template.simple.bed')
        bed = pd.read_csv(self.simple_template,header=None, sep='\t')
        chr_list = ['chr{0}'.format(a) for a in list(range(1,23))+['X','Y']]
        bed = bed[bed[0].isin(chr_list)]
        chr_weight = bed[0].value_counts()/len(bed)
        len_weight = bed[3].value_counts()/len(bed)
        return chr_weight.to_frame(name='weight'), len_weight.to_frame(name='weight')
    
    def cal_chimeric_weight(self):
        '''
        calculate the chr and length weight of compelx eccDNA  
        '''
        if not self.chimeric_template:
            self.chimeric_template = join(split(__file__)[0], 'resource', 'template', 'template.chimeric.bed')
            ## calculate from template
        bed = pd.read_csv(self.chimeric_template,header=None, sep='\t')
        fragN_sum = bed[4].value_counts()
        fragN_weight = fragN_sum[fragN_sum<10].value_counts()/len(fragN_sum[fragN_sum<10])
        chr_weight = {}; len_weight  = {}
        for i in range(2,10):
            temp = bed[bed[4].isin(list(fragN_sum[fragN_sum==i].index))]
            chr_weight[i] = temp[0].value_counts().to_frame(name='weight') / temp[0].value_counts().sum()
            len_weight[i] = temp[3].value_counts().to_frame(name='weight') / temp[3].value_counts().sum()
        return fragN_weight.to_frame(name='weight'), chr_weight, len_weight
    
class seqsim(object):
    '''
    sample: sample name
    reference: path of reference fasta 
    circular_number: total number of circular DNA to simulate
    linear_number: total number of circular DNA to simulate
    seed: srandom seed to reproduce the result
    path: path to preserve result files
    simple_ratio: the ratio of simple eccDNA in simulated data (optional)
    simple_template: bed file of simple eccDNA in format 'chr start end length id'
    chimeric_template: bed file of chimeric eccDNA in format 'chr start end length id'
    '''
    def __init__(self, 
                 sample: str,
                 reference: str,
                 path: str,
                 linear_number: int,
                 circular_number: int,
                 seed = None,
                 simple_ratio = None,
                 simple_template = None, 
                 chimeric_template = None,
                ):
        
        self.sample = [sample, [sample+'_sn', sample+'_sp'], [sample+'_cn', sample+'_cp'], [sample+'.neg', sample+'.pos']]
        self.number = [linear_number, circular_number]
        self.seed = seed
        self.path = join(path, self.sample[0])
        self.utils = utilities(reference)
        self.prof = statistics(simple_ratio, simple_template, chimeric_template)
        self.prefix = [[join(self.path, self.sample[a][0]), join(self.path, self.sample[a][1])] for a in [1,2,3]]    
        
        if not exists(self.path):
            os.makedirs(self.path)
        self.sim_sequence()

    def random_simple(self):
        '''
        function to simulate simple DNA
        '''
        np.random.seed(self.seed)
        
        ## calculate the weight matrix and genome length table
        weight = self.prof.simple_weight
        reference_len = self.utils.genome_length()

        ##random generate eccDNA
        for type_i in range(2):
            ## init the result dataframe
            output_table = pd.DataFrame(columns=['id', 'fragN','region', 'length', 'seq'])
            output_bed = pd.DataFrame(columns=['chrom', 'start', 'end', 'length', 'id'])
            for ecc_i in range(int(self.number[type_i]*self.prof.simple_ratio)):
                chrom = np.random.choice(weight[0].index, p = weight[0]['weight'])
                length = np.random.choice(weight[1].index, p = weight[1]['weight'])
                while (length > reference_len.loc[chrom]['length']):
                    length = np.random.choice(weight[1].index, p = weight[1]['weight'])
                ecc = self.utils.random_region(chrom, length)
                ## transfer into output format
                region = ecc[0] + ':' + str(ecc[1]) + '-' + str(ecc[2])
                output_bed.loc[ecc_i] = ecc[:4] + [self.sample[1][type_i] + str(ecc_i+1)]
                output_table.loc[ecc_i] = [self.sample[1][type_i] + str(ecc_i+1)] + [1 , region, ecc[3], ecc[4]]    
            ## write result file
            self.utils.write_fasta(output_table, self.prefix[0][type_i] + '.fasta')
            output_table.to_csv(self.prefix[0][type_i] + '.csv', index=None, sep='\t')
            output_bed.to_csv(self.prefix[0][type_i] + '.bed', header=None, index=None, sep='\t')
        return
    
       
    def random_chimeric(self):
        '''
        function to simulate chimeric DNA
        '''
        np.random.seed(self.seed)
        ## calculate the weight matrix and genome length table
        weight = self.prof.chimeric_weight
        reference_len = self.utils.genome_length()
        
        for type_i in range(2):
            ## init result dataframe
            output_table = pd.DataFrame(columns= ['id', 'fragN', 'region', 'length', 'seq'])
            output_bed = pd.DataFrame(columns= ['chrom', 'start', 'end', 'length', 'id'])
            frag_i = 0
            for ecc_interger in range(int(self.number[type_i]*self.prof.simple_ratio), self.number[type_i]):
                ## initiate ecc
                fragN = np.random.choice(weight[0].index, p = weight[0]['weight'])
                ecc_seq = ''; ecc_region = []; ecc_length = 0
                ## random generate fragment
                for i in range(fragN):
                    frag_chrom = np.random.choice(weight[1][fragN].index, p = weight[1][fragN]['weight'])
                    frag_length = np.random.choice(weight[2][fragN].index, p = weight[2][fragN]['weight'])
                    while (frag_length > reference_len.loc[frag_chrom]['length']):
                        frag_length = np.random.choice(weight[1].index, p = weight[1]['weight'])
                    frag = self.utils.random_region(frag_chrom, frag_length); 
                    frag_i += 1
                    ## transfer into output format
                    output_bed.loc[frag_i] = frag[0:4] + [''.join([self.sample[2][type_i], str(ecc_interger+1)])]
                    ecc_seq += frag[4]; ecc_region += [frag[0] + ':' + str(frag[1]) + '-' + str(frag[2])] ; ecc_length += frag[3]
                ecc_region = '|'.join(ecc_region)
                output_table.loc[ecc_interger] = [''.join([self.sample[2][type_i], str(ecc_interger+1)])] + [fragN, ecc_region, ecc_length, ecc_seq]
                
            ## write result file
            self.utils.write_fasta(output_table, self.prefix[1][type_i] + '.fasta')
            output_table.to_csv(self.prefix[1][type_i] + '.csv', index=None, sep='\t')
            output_bed.to_csv(self.prefix[1][type_i] + '.bed', header=None, index=None, sep='\t')
        return 
        
    def merge(self):
        '''
        function to merge simulated simple & complex file
        '''
        for type_i in range(2):
            if (self.prof.simple_ratio==1):
                simple_table = pd.read_csv(self.prefix[0][type_i] + '.csv', sep='\t')
                simple_table.to_csv(self.prefix[2][type_i] + '.csv', index=None, sep='\t')
                for format_ in ['.bed', '.fasta']:
                    merge_cmd = ' '.join(['cat', 
                                          self.prefix[0][type_i] + format_,
                                          '>', 
                                          self.prefix[2][type_i] + format_])
                    sp.check_call(merge_cmd, shell=True)
                clean_cmd = 'rm {0}*'.format(self.prefix[0][type_i])
                sp.check_call(clean_cmd, shell=True)
            if (self.prof.simple_ratio==0):
                chimeric_table = pd.read_csv(self.prefix[1][type_i] + '.csv', sep='\t')
                chimeric_table.to_csv(self.prefix[2][type_i] + '.csv', index=None, sep='\t')
                for format_ in ['.bed', '.fasta']:
                    merge_cmd = ' '.join(['cat', 
                                          self.prefix[1][type_i] + format_,
                                          '>', 
                                          self.prefix[2][type_i] + format_])
                    sp.check_call(merge_cmd, shell=True)
                clean_cmd = 'rm {0}*'.format(self.prefix[1][type_i])
                sp.check_call(clean_cmd, shell=True)
            if (self.prof.simple_ratio!=0)&(self.prof.simple_ratio!=1):
                simple_table = pd.read_csv(self.prefix[0][type_i] + '.csv', sep='\t')
                chimeric_table = pd.read_csv(self.prefix[1][type_i] + '.csv', sep='\t')
                pd.concat([simple_table, chimeric_table]).to_csv(self.prefix[2][type_i] + '.csv', index=None, sep='\t')
                for format_ in ['.bed', '.fasta']:
                    merge_cmd = ' '.join(['cat', 
                                          self.prefix[0][type_i] + format_,
                                          self.prefix[1][type_i] + format_,
                                          '>', 
                                          self.prefix[2][type_i] + format_])
                    sp.check_call(merge_cmd, shell=True)
                clean_cmd = 'rm {0}* {1}*'.format(self.prefix[0][type_i],self.prefix[1][type_i])
                sp.check_call(clean_cmd, shell=True)
        return
    
    def sim_sequence(self):
        if (self.prof.simple_ratio!=0):
            self.random_simple()
        if (self.prof.simple_ratio!=1):
            self.random_chimeric()
        self.merge()
        return
    
class libsim():
    '''
    sample: sample name
    reference: path of reference fasta 
    seed: srandom seed to reproduce the result
    amp: length(bp) of RCA
    path: path to preserve result files
    '''
    def __init__(self, 
                 sample: str, 
                 reference : str,
                 path: str, 
                 seed = None,
                 meancov = 25,
                 amp = 5000, 
                ):
        
        self.sample = sample
        
        self.seed = seed
        self.meancov = meancov
        self.amp = amp

        self.path = join(path, self.sample)
        self.prefix = [join(self.path, self.sample + '.neg'), join(self.path, self.sample + '.pos'), join(self.path,  self.sample + '.lib')]
        self.utils = utilities(reference)
        
        self.sim_library()
             
    def sim_breakpoint(self, region_table):
        '''
        random generate shift to simulate breakpoint and return the shifted region_table
        '''
        np.random.seed(self.seed)
        region_table = region_table.reset_index().drop('index', axis=1)
        shifted_table = pd.DataFrame(columns=['chrom', 'start', 'end', 'length', 'seq', 'id'])
        shift = np.random.randint(1, region_table.length.sum())
        for i in range(1, len(region_table) + 1): ## residual
            if (shift == region_table[0:i].length.sum()):
                shifted_table = pd.concat([region_table[i:].copy(),region_table[0:i].copy()],axis=0).reset_index().drop('index', axis=1)
                break
            elif (shift < region_table[0:i].length.sum()):
                shift_frag = shift - region_table[0:i-1].length.sum()
                shifted_table = pd.concat([region_table[i-1:].copy(),region_table[0:i].copy()],axis=0).reset_index().drop('index', axis=1)
                shifted_table.loc[0, 'start'] += shift_frag
                shifted_table.loc[0, 'length'] -= shift_frag
                shifted_table.loc[0, 'seq'] = region_table.loc[i-1, 'seq'][shift_frag:]
                shifted_table.loc[len(shifted_table)-1, 'end'] = region_table.loc[i-1, 'start'] + shift_frag
                shifted_table.loc[len(shifted_table)-1, 'length'] = shift - region_table[0:i-1].length.sum()
                shifted_table.loc[len(shifted_table)-1, 'seq'] = region_table.loc[i-1, 'seq'][:shift_frag]
                break
        return shifted_table
    
    def sim_RCA(self, region_table):
        '''
        amplify to simulate RCA
        '''
        ## initiate the amplified table
        region_table = region_table.reset_index().drop('index', axis=1)
        amplified_table = pd.DataFrame(columns=['chrom', 'start', 'end', 'length', 'seq', 'id'])
        ## amplified 
        roundN = int(self.amp/ region_table.length.sum())
        for i in range(roundN + 1): ## round
            amplified_table = pd.concat([amplified_table,region_table],axis=0)
        residual = self.amp - roundN * region_table.length.sum()
        for i in range(len(region_table)+1): ## residual
            if region_table[0:i].length.sum()>residual:
                residual_table = region_table[0:i].copy()
                residual_table.loc[i-1,'length'] = residual - residual_table[0:i-1].length.sum()
                residual_table.loc[i-1,'end'] = int(residual_table.loc[i-1,'start']) + int(residual_table.loc[i-1,'length'])
                residual_table.loc[i-1,'seq'] = residual_table.loc[i-1,'seq'][0:(residual - residual_table[0:i-1].length.sum())]
                break
        amplified_table = pd.concat([amplified_table, residual_table],axis=0)
        amplified_table = amplified_table.reset_index()
        amplified_table = amplified_table.drop('index',axis=1)
        return amplified_table
    
    def sim_library(self):
        '''
        '''
        output_table = pd.DataFrame(columns= ['id', 'fragN', 'region', 'length', 'seq'])
        output_bed = pd.DataFrame(columns= ['chrom', 'start', 'end', 'length', 'id'])
        ## positive library
        pos_bed = pd.read_csv(self.prefix[1] + '.bed', sep='\t',names=['chrom', 'start', 'end', 'length', 'id'])
        pos_bed.insert(pos_bed.shape[1], 'seq', pos_bed.apply(lambda x: self.utils.get_seq(x.chrom, x.start, x.end), axis=1))
        if (self.amp==0):
            for i in pos_bed.id.unique():
                processed_ecc = self.sim_breakpoint(pos_bed[pos_bed.id == i])
                processed_ecc['region'] = processed_ecc['chrom'] + ':' + processed_ecc['start'].astype(str) + '-' + processed_ecc['end'].astype(str)
                output_table.loc[i] = [i, len(processed_ecc), '|'.join(processed_ecc['region']), processed_ecc['length'].sum(), ''.join(processed_ecc['seq'])]
                output_bed= pd.concat([output_bed, processed_ecc[['chrom', 'start', 'end', 'length', 'id']]])
        else:
            for i in pos_bed.id.unique():
                processed_ecc = self.sim_breakpoint(pos_bed[pos_bed.id == i])
                processed_ecc = self.sim_RCA(processed_ecc)
                processed_ecc['region'] = processed_ecc['chrom'] + ':' + processed_ecc['start'].astype(str) + '-' + processed_ecc['end'].astype(str)
                output_table.loc[i] = [i, len(processed_ecc), '|'.join(processed_ecc['region']), processed_ecc['length'].sum(), ''.join(processed_ecc['seq'])]
                output_bed= pd.concat([output_bed, processed_ecc[['chrom', 'start', 'end', 'length', 'id']]])
        ## merge negative library & positive library
        neg_table = pd.read_csv(self.prefix[0] + '.csv', sep='\t')
        output_table = pd.concat([output_table, neg_table])
        output_table.insert(output_table.shape[1],'realcov',np.random.gamma(2.5, self.meancov/2.5, len(output_table)))
        output_table.insert(output_table.shape[1],'tempcov',output_table.realcov*(output_table.length)/(output_table.length+self.amp))
        output_table.to_csv(self.prefix[2] + '.csv', index=None, sep='\t')
        self.utils.write_fasta(output_table, self.prefix[2] + '.fasta')
        output_bed.to_csv(self.prefix[2] + '.bed', header=None, index=None, sep='\t')
        for format_ in ['.bed', '.fasta']:
            merge_cmd = ' '.join(['cat', self.prefix[0] + format_, '>>', self.prefix[2] + format_])
            sp.check_call(merge_cmd, shell=True)
        return
    

class fqsim(object):
    '''
    sample: sample name
    csv: path of template csv file generated in lib sim
    path: path of result
    
    thread: maximum thread to run 
    meancov: mean coverage of eccDNA
    ont_mean
    ont_std
    ont_model
    sr_mean
    sr_std
    sr_readlen
    sr_platform
    '''
    def __init__(self, 
                 sample: str,
                 csv: str,
                 path: str,
                 thread = 8,
                 ont_mean = 3000,
                 ont_std = 2500,
                 ont_model = 'R94',
                 sr_mean = 400,
                 sr_std = 125,
                 sr_readlen = 150,
                 sr_platform = 'HS25'
                ):
        
        self.sample = sample
        self.path = join(path, self.sample)
        self.thread = int(thread)
        self.ont_mean = float(ont_mean)
        self.ont_std = float(ont_std)
        self.ont_model = ont_model
        self.sr_mean = float(sr_mean)
        self.sr_std = float(sr_std)
        self.sr_readlen = int(sr_readlen)
        self.sr_platform = sr_platform
        
        self.utils = utilities('')
        self.ec = pd.read_csv(csv,sep='\t')
        self.unifa = join(self.path, 'unifa')
        self.tmp = join(self.path, 'tmp')
        self.tool_path = split(__file__)[0]
        
        if not exists(self.unifa):
            os.makedirs(self.unifa) 
        if not exists(self.tmp):
            os.makedirs(self.tmp) 
        
        self.unifasta()
        self.sim_fastq_sr()
        self.sim_fastq_ont()
        os.chdir(self.path)
        self.sort()
        sp.check_call('rm -rf {0} {1}'.format(self.unifa, self.tmp), shell=True)
        
            
    def unifasta(self):
        '''
        Function to write each single eccDNA into a fasta file
        '''
        for i in self.ec.index:
            tmp_fasta = join(self.unifa, self.ec.loc[i,'id'] + '.fasta')
            self.utils.write_fasta(self.ec.loc[i:i,:], tmp_fasta)
        return
    
    def multi_para_art(self):
        para_file = [] 
        for i in self.ec.index:
            input_ = join(self.unifa, self.ec.loc[i,'id'] + '.fasta')
            output_ = join(self.tmp, self.ec.loc[i,'id'] + '.R')
            para_file.append((input_, output_, self.sr_platform, self.sr_readlen, self.ec.loc[i,'tempcov'], self.sr_mean, self.sr_std))
        return para_file
    
    def art(self, input_, output_, sr_platform, length, cov, fral, frastd):
        art_cmd = 'art_illumina -na -q -nf 0 -p -i {0} -o {1} -ss {2} -l {3} -f {4} -m {5} -s {6}'.format(input_, output_, sr_platform, length, cov, fral, frastd)
        sp.check_call(art_cmd, shell=True)
        self.utils.transfer_files('{0}1.fq'.format(output_), '{0}/{1}.sr.R1.fastq'.format(self.path, self.sample))
        self.utils.transfer_files('{0}2.fq'.format(output_), '{0}/{1}.sr.R2.fastq'.format(self.path, self.sample))
        clean_cmd = 'rm {0}*'.format(output_)
        sp.check_call(clean_cmd, shell=True)
        return
    
    def sim_fastq_sr(self):
        '''
        funtions to generate NGS fastq reads
        '''
        self.art_para_file = self.multi_para_art()
        os.chdir(self.tmp)
        p = mp.Pool(self.thread)
        for i in range(int(len(self.art_para_file)/self.thread)+1):
            p.starmap_async(self.art, self.art_para_file[i*self.thread:(i+1)*self.thread])
        p.close()
        p.join()
        return
    
    def multi_para_pbsim2(self):
        para_file = []
        model_file = join(self.tool_path, 'resource', 'pbsim2', self.ont_model + '.model')
        for i in self.ec.index:
            input_ = join(self.unifa, self.ec.loc[i,'id'] + '.fasta')
            output_ = join(self.tmp, self.ec.loc[i,'id'])
            para_file.append((self.ec.loc[i,'tempcov'], output_, self.ec.loc[i,'id'], model_file, self.ont_mean, self.ont_std, input_))
        return para_file
    
    def pbsim2(self, cov, output_, id_prefix, model, ont_mean, ont_std, input_):
        cmd = 'pbsim  --depth {0} --prefix {1} --id-prefix {2} --hmm_model {3} --length-mean {4} --length-sd {5} {6}'.format(cov, output_, id_prefix, model, ont_mean, ont_std, input_)
        sp.check_call(cmd, shell=True)
        self.utils.transfer_files('{0}_0001.fastq'.format(output_), '{0}/{1}.ont.fastq'.format(self.path, self.sample))
        clean_cmd = 'rm {0}*'.format(output_)
        sp.check_call(clean_cmd, shell=True)
        return
    
    def sim_fastq_ont(self):
        '''
        funtions to generate ONT fastq reads
        '''
        os.chdir(self.tmp)
        p = mp.Pool(self.thread)
        self.pbsim2_para_file = self.multi_para_pbsim2()
        for i in range(int((len(self.pbsim2_para_file)-1)/self.thread)+1):
            p.starmap_async(self.pbsim2, self.pbsim2_para_file[(i*self.thread+1):(i+1)*self.thread+1])
        p.close()
        p.join()
        return
    
    def sort(self):
        '''
        functions to sort fastq
        '''
        sort_cmd1 = 'seqkit sort {0}/{1}.sr.R1.fastq > {0}/{1}.sr.sorted.R1.fastq '.format(self.path, self.sample)
        sort_cmd2 = 'seqkit sort {0}/{1}.sr.R2.fastq > {0}/{1}.sr.sorted.R2.fastq '.format(self.path, self.sample)
        sort_cmd3 = 'seqkit sort {0}/{1}.ont.fastq > {0}/{1}.ont.sorted.fastq '.format(self.path, self.sample)
        clean_cmd = 'rm -rf {0}/{1}.ont.fastq {0}/{1}.sr.R1.fastq  {0}/{1}.sr.R2.fastq {0}/unifa {0}/tmp'.format(self.path, self.sample)
        sp.check_call(sort_cmd1, shell=True)
        sp.check_call(sort_cmd2, shell=True)
        sp.check_call(sort_cmd3, shell=True)
        sp.check_call(clean_cmd, shell=True)
        return
    
