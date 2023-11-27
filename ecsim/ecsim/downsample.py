
import os
import sys
import argparse
from os.path import exists,join
from importlib.resources import files
from .simulate import seqsim, libsim, fqsim
from .utils import utilities

class sim_all(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser(prog='ecsim',
                                 description='generating simulated dataset',
                                 epilog='Text at the bottom of help'
                                )
        ### Required Parameters
        self.parser.add_argument('--sample', required=True, type=str, help='sample name')
        self.parser.add_argument('--reference', required=True, type=str, help='path of reference fasta file')
        self.parser.add_argument('--thread', required=True, type=int, help='maximum thread to run')
        self.parser.add_argument('--path', required=True, type=str, help='output path of simulated data')
        
        ## Optional Parameters
        self.parser.add_argument('--meancov', required=False, type=float, default=25, help='mean coverage of all DNA molecules. (default=25)')
        self.parser.add_argument('--circular-number', required=False, type=int, default=5000, help='number of circular DNA molecules to simulate. (default=5000)')
        self.parser.add_argument('--linear-number', required=False, type=int, default=5000, help='number of linear DNA molecules to simulate. (default=5000)')
        self.parser.add_argument('--amp', required=False, type=int, default=5000, help='length(bp) of RCA. (default=5000)')
        self.parser.add_argument('--seed', required=False, type=int, default=None, help='set a random seed to repeat result.(default=None)')
        
        ### Simple & chimeric DNA sample file
        self.parser.add_argument('--simple-ratio', required=False, type=float, help='ratio of simple DNA to simulate. (default is calculate from the profile)')
        self.parser.add_argument('--simple-profile', required=False, type=str, help='simple profile (format: chr start end length id) to simulate from. (alternative)')
        self.parser.add_argument('--chimeric-profile', required=False, type=str, help='chimeric profile (format: chr start end length id) to simulate from. (alternative)')
        
        ### Fastq parameters
        self.parser.add_argument('--ont-model', required=False, type=str, default='R94', help='paras for PBSIM2 option(R94, R95, R103, P4C2, P5C3, P6C4)')
        self.parser.add_argument('--ont-mean', required=False, type=float, default=3000, help='mean of reads length of simulated ONT fastq (default=3000)')
        self.parser.add_argument('--ont-std', required=False, type=float, default=2500, help='std  of reads length of simulated ONT fastq (default=2500)')
        self.parser.add_argument('--sr-platform', required=False, type=str, default='HS25', help='paras for art option(HS10, HS20,HS25, HSXn, HSXt, MinS, MSv1, MSv3, NS50)')
        self.parser.add_argument('--sr-mean', required=False, type=float, default=400, help='mean of insert length of simulated NGS fastq (default=400)')
        self.parser.add_argument('--sr-std', required=False, type=float, default=125, help='std  of insert length of simulated NGS fastq (default=125)')
        self.parser.add_argument('--sr-readlen', required=False, type=float, default=150, help='reads length of simulated NGS fastq (default=125)')
        self.args = self.parser.parse_args()
        
        self.simulate()
            
    def simulate(self): 
        temp = fqsim(sample = self.args.sample,
                     csv = join(self.args.path, self.args.sample, self.args.sample + '.lib.csv'),
                     path = self.args.path,
                     thread = self.args.thread,
                     ont_model = self.args.ont_model,
                     ont_mean = self.args.ont_mean,
                     ont_std = self.args.ont_std,
                     sr_platform = self.args.sr_platform,
                     sr_mean = self.args.sr_mean,
                     sr_std = self.args.sr_std,
                     sr_readlen = self.args.sr_readlen,
                     )

def main():
    run = sim_all()
    
if __name__ == '__main__':
    main()
