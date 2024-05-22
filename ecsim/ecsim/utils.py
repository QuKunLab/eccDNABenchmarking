#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import pysam as ps 
import subprocess as sp


class utilities(object):
    def __init__(self, 
                 reference: str,
                ):
        """
        utilities
        """
        self.reference = reference
        
    def genome_length(self):
        """
        calculate the length of each chromosome
        """
        fasta = ps.Fastafile(self.reference)
        output = pd.DataFrame([fasta.lengths], columns=fasta.references, index=['length']).T
        output = output.drop(list(output.filter(regex='_',axis=0).index)+['chrM'])
        return output
    
    def get_seq(self, chrom, start, end):
        fasta = ps.Fastafile(self.reference)
        seq = fasta.fetch(chrom, start, end)
        seq = seq.upper()
        return seq
    
    def random_region(self, chrom, length):
        genome = self.genome_length()
        start = np.random.randint(genome.loc[chrom]['length']-length)
        end = start + length
        seq = self.get_seq(chrom, start, end)
        output = [chrom, start, end, length, seq]
        if (seq.count('N')>0):
            output = self.random_region(chrom, length)
        return output
    
    def transfer_files(self, file1, file2):
        f1 = open(file1,'r').read()        
        with open(file2,'a') as f2:
            f2.write(f1)
        f2.close()
        return
    
    def write_fasta(self, data, output):
        with open(output,'w') as f:
            for i in data.index:
                f.write('>{0}\n'.format(data.loc[i,'id']))
                f.write(data.loc[i,'seq'])
                f.write('\n')
            f.close()
        return
    
