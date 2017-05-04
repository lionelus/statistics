from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import bz2
import csv
import re
import glob
from matplotlib import pyplot as plt
from collections import Counter
import pandas as pd

def get_MinknowVersion(h5py_file):
    """
    Get the Minknow version from fast5 file
    """
    version = list(h5py_file['/UniqueGlobalKey/tracking_id'].attrs.items())
    version_d = {key: value.decode('utf-8') for key, value in version}
    return version_d['version']
#print(get_MinknowVersion())

def getFlowcellId(h5py_file):
    """
    Get the flowcell id from fast5 file
    """
    flowcell_id = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    flowcell_id_dico = {key: value.decode('utf-8') for key, value in flowcell_id}
    return flowcell_id_dico['flow_cell_id']

def get_Hostname(h5py_file):
    """
    Get the hostname from fast5 file
    """
    host_name = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    host_name_dico = {key: value.decode('utf-8') for key, value in host_name}
    return host_name_dico['hostname']

def getNumMinION(h5py_file):
    """
    Get the number of Minion run
    """
    numMinION = list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    numMinION_dico = {key: value.decode('utf-8') for key, value in numMinION}
    return numMinION_dico['device_id']

def getProtocolRunId(h5py_file):
    """
    Get the run id protocol from fast 5 file
    """
    protocol_run_id =  list(h5py_file["/UniqueGlobalKey/tracking_id"].attrs.items())
    protocol_run_id_dico = {key: value.decode('utf-8') for key, value in protocol_run_id}
    return protocol_run_id_dico['protocol_run_id']

def get_barcode():
    """
    Get the barcode from a file given in input
    """
    barcode_file = input('Enter the name of the file where the barcode is stored:')
    set_doublon = set()

    with open(barcode_file) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')

        for row in spamreader:
            pattern = re.search(r'BC(\d{2})', row[0])

            if pattern:
                barcode = 'barcode{}'.format(pattern.group(1))
                set_doublon.add(barcode)
    return list(set_doublon)


#selection : barcode name list
def get_fastq(selection):
    """
    Get the fastq sequence
    """
    global_length_array = []
    path_bz2_file = '/home/ferrato/Documents/salut.txt'
    path_bz2_directory = input('path to bz2 directory:')

    if os.path.isfile(path_bz2_file):
        open(path_bz2_file, 'w').close()

    if not os.path.exists('images'):
        os.makedirs('images')

    if os.path.exists('images'):
        for file in os.listdir('images'):
            os.remove('images/'+file)

    if not os.path.exists('statistics'):
        os.makedirs('statistics')

    if os.path.exists('statistics'):
        for file in os.listdir('statistics'):
            os.remove('statistics/'+file)

    for element in glob.glob("{}/*.bz2".format(path_bz2_directory)):

        counter_template = Counter()
        total_nucs_template = 0
        barcode_length_array = []

        for el in selection:
            if el in element:
                print(element)
                uncompressedData = bz2.BZ2File(element,'rb').read()
                uncomp = uncompressedData.decode('utf-8')
                file = open(path_bz2_file,'w')
                file.write(uncomp)

                with open(path_bz2_file) as in_handle:

                    for title, seq, qual in FastqGeneralIterator(in_handle):
                        barcode_length_array.append(len(seq))
                        global_length_array.append(len(seq))

                        for nucleotide_template in seq:
                            counter_template[nucleotide_template] += 1
                            total_nucs_template += 1

                completeName = os.path.join('statistics/',el)
                barcode_file = open(completeName,'w')

                series_read_size = pd.Series(barcode_length_array)
                statistics = pd.Series.describe(series_read_size)

                for index, value in statistics.iteritems():
                    barcode_file.write("Read.fastq.length.{}={}\n".format(index, value))

                for nucleotide, count in counter_template.items():
                    barcode_file.write("nucleotide.{}.template={}\n".format(nucleotide, count))
                    if nucleotide == 'total':
                        continue
                    calcul = float(count) / float(total_nucs_template)
                    barcode_file.write("nucleotide.{}.proportion={}\n".format(nucleotide, calcul))
                barcode_file.close()

                plt.boxplot(barcode_length_array, showfliers=False)
                plt.title('Boxplot of read length by barcode')
                plt.savefig('images/image{}.png'.format(el))
                plt.close()
                file.close()


#    with open(path_bz2_file) as in_handle:
#
 #       for title, seq, qual in FastqGeneralIterator(in_handle):
  #          length_array.append(len(seq))
#
 #           for nucleotide_template in seq:
  #              counter_template[nucleotide_template] += 1
  #              total_nucs_template += 1


  #  if length_array == []:
  #      print('Le dossier indiqué ne contient pas de fichiers bz2')


    return global_length_array



