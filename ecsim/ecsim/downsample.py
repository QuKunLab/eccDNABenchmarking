import os
import sys
import argparse
from os.path import exists,join
from .simulate import seqsim, libsim, fqsim
from .utils import utilities

class sim_all(object):

    def __init__(self):
        self.parser = argparse.ArgumentParser(prog='downsample',
                                 description='generating fastq using the downsampled csv',
                                 epilog='_______________________________________________'
                                )
        ### Required Parameters
        self.parser.add_argument('--sample', required=True, type=str, help='sample name (director containing downsampled csv in the path)')
        self.parser.add_argument('--thread', required=True, type=int, help='maximum thread to run')
        self.parser.add_argument('--path', required=True, type=str, help='output path of simulated data')

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
