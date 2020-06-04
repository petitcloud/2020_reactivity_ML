# Calculating Fukui distribution factor (FDF) from parr.txt and index_list.txt
import sys, string, pprint
import numpy
import pandas as pd
from math import log10

# Read index_list.txt

def gen_parr_array(dirc):

 indexfile = open('{0}/total_index_list.txt' .format(dirc), 'r')
 indexline = indexfile.readlines()

 folder_name = []
 index_value = []
# exp_elec = []

 for a in range(len(indexline)):
  folder_name.append(indexline[a].split()[0])
  index_value.append(float(indexline[a].split()[2]))
#  exp_elec.append(float(indexline[a].split()[2]))

# exp_elec = numpy.array(exp_elec)

 # Read parr.txt in each folder and calculate FDF

 parr_array = []
 parr_ind_array = []
 len_array =[]

 for b in range(len(folder_name)):
  tempfile = pd.read_excel('{0}/{1}/ea_adia/Parr.xlsx' .format(dirc, folder_name[b]), 'Sheet1')

#  tempfile = open('{0}/{1}/parr.txt' .format(dirc, folder_name[b]), 'r')
#  templines = tempfile.readlines()
#  num_atom = int((len(templines)-4)/2)
#  temp_parr = []
#  temp_array = []

#  for bb in range(num_atom):
#   temp_parr.append(float(templines[bb-num_atom].split()[2]))
 
  temp_array = numpy.array(tempfile['Spin density'])
  parr_array.append(temp_array)
  parr_ind_array.append(list(index_value[b]*temp_array)) 
  len_array.append(len(tempfile['Spin density']))

 parr_mat = numpy.zeros(len(folder_name)*max(len_array)).reshape(len(folder_name),max(len_array))
 parr_ind_mat = numpy.zeros(len(folder_name)*max(len_array)).reshape(len(folder_name),max(len_array))
 for c in range(len(folder_name)):
  for cc in range(len_array[c]):
   parr_mat[c,cc] = pow(10, parr_array[c][cc])
   parr_ind_mat[c,cc] = pow(10, parr_ind_array[c][cc])

 return folder_name, parr_mat, parr_ind_mat
 
# Print results
if __name__ == '__main__':
 folder_name, parr_mat, parr_ind_mat = gen_parr_array('/home/petitcloud/electrolyte_air_screen')
 print ('<folder_name>\t<ea*p1>\t<pdf>\t')
 for d in range(len(folder_name)):
  print (folder_name[d], log10(parr_ind_mat[d][0]), log10(sum(parr_ind_mat[d])))
#  print(folder_name[d], parr_mat[d])

