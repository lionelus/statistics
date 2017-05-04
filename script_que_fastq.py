from matplotlib import pyplot as plt
import pandas as pd
from collections import Counter
import os
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import getter1D
import numpy as np
import bz2

class fastq:
    def __init__(self):
        #self.selection = getter1D.get_barcode()
        #self.selection.sort()
        self.fastq_length_array, self.counter_template, self.total_nucs_template = get_fastq_bz2()

    def statistics_read_size(self):
        """
        Get statistics from the file containing the fastq bz2 files decompressed
        """
        series_read_size = pd.Series(self.fastq_length_array)
        statistic1 = pd.Series.describe(series_read_size)
        return statistic1

    def reads_size(self):

        if self.fastq_length_array == []:
            print('There is a mistake')
            return False
        plt.hist(self.fastq_length_array, edgecolor="#E6E6E6", color="#EE6666", bins=range(min(self.fastq_length_array), max(self.fastq_length_array) + 100, 100))
        plt.xlabel("fastq size")
        plt.ylabel("Count")
        plt.xlim(0,6000)
        plt.title("reads size")
        plt.savefig('images/image7.png')
        plt.close()

    def counter(self):
        """
        Paricipate to the count df nucleotide(A,T,C,G)" \
        """
        return self.counter_template, self.total_nucs_template

def log_file1D(basecall_stat):

    counter_template, total_nucleotide_template = basecall_stat.counter()

    completeName = os.path.join('/home/ferrato/Documents/fast5', "fichier_aozan.txt")

    with open(completeName, 'w') as file_data:

        for nucleotide, count in counter_template.items():
            file_data.write("nucleotide.{}.template={}\n".format(nucleotide, count))


        for index, element in basecall_stat.statistics_read_size().iteritems():
            file_data.write("Read.fastq.length.{}={}\n".format(index, element))


def get_fastq_bz2():
    """
    Get the fastq sequence
    """
    counter_template = Counter()
    total_nucs_template = 0

    length_array = []


    path_bz2_directory = input('path to bz2 file:')


    with open(path_bz2_directory) as in_handle:

        for title, seq, qual in FastqGeneralIterator(in_handle):
            length_array.append(len(seq))

            for nucleotide_template in seq:
                counter_template[nucleotide_template] += 1
                total_nucs_template += 1


        if length_array == []:
            print('Le dossier indiqué ne contient pas de fichiers bz2')

    return length_array, counter_template, total_nucs_template











  #  for element in glob.glob("{}/*.bz2".format(path_bz2_directory)):

  #      for el in selection:
   #         if el in element:
    #            print(element)
     #           uncompressedData = bz2.BZ2File(element, 'rb').read()
      #          uncomp = uncompressedData.decode('utf-8')
       #         file = open(path_bz2_file, 'a')
        #        file.write(uncomp)
         #       file.close()



    #with open(path_bz2_file) as in_handle:

     #   for title, seq, qual in FastqGeneralIterator(in_handle):
      #      length_array.append(len(seq))

       #     for nucleotide_template in seq:
        #        counter_template[nucleotide_template] += 1
         #       total_nucs_template += 1


      #  if length_array == []:
      #      print('Le dossier indiqué ne contient pas de fichiers bz2')

    #return length_array, counter_template, total_nucs_template

fast = fastq()
log_file1D(fast)

#path_bz2_file = '/home/ferrato/Documents/salut.txt'
#path_bz2_directory = input('path to bz2 file:')

#if os.path.isfile(path_bz2_file):
 #   open(path_bz2_file, 'w').close()

#uncompressedData = bz2.BZ2File(path_bz2_directory, 'rb').read()
#uncomp = uncompressedData.decode('utf-8')
#file = open(path_bz2_file, 'a')
#file.write(uncomp)
#file.close()