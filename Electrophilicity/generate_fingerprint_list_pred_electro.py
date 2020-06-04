
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem import AllChem

from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import Molecule

import numpy as np
import os


# In[2]:


class gen_fingerprint_array:

    FILE_NAME = 'redox'
    folder_name = []
#    exp_nuc = []
    smiles = []
    fp_MAC = []
    fp_TOP = []
    fp_MOR = []   

    def __init__(self, dirc):
        with open('{}/total_index_list.txt' .format(dirc), 'r') as indexfile:
            indexline = indexfile.readlines()

        for a in range(len(indexline)):
            self.folder_name.append(indexline[a].split()[0])
#            self.exp_nuc.append(float(indexline[a].split()[2]))
#            self.smiles.append(indexline[a].split()[3])
        
#        self.exp_nuc = np.array(self.exp_nuc)
        
        for b in range(len(indexline)):
            if not os.path.exists('{}/{}/temp.mol' .format(dirc, self.folder_name[b])):
                try:
                    try:
                        tempfile = GaussianOutput('{}/{}/{}.log' .format(dirc, self.folder_name[b], self.FILE_NAME))
                    except:
                        tempfile = GaussianOutput('{}/{}/ea_adia/{}.log' .format(dirc, self.folder_name[b], self.FILE_NAME))
                except:
                    tempfile = GaussianOutput('{}/{}/vertical_ie/{}.log' .format(dirc, self.folder_name[b], self.FILE_NAME))
                tempmol = tempfile.final_structure
                tempmol.to(fmt='mol', filename='{}/{}/temp.mol' .format(dirc, self.folder_name[b]))
            rdtemp = Chem.MolFromMolFile('{}/{}/temp.mol' .format(dirc, self.folder_name[b]))
#            if rdtemp is None:
#                rdtemp = Chem.MolFromSmiles(self.smiles[b])
            rd_mac = MACCSkeys.GenMACCSKeys(rdtemp)
            rd_top = RDKFingerprint(rdtemp)
            rd_mor = AllChem.GetMorganFingerprintAsBitVect(rdtemp, 2)
            
            self.fp_MAC.append(np.array([int(x) for x in rd_mac.ToBitString()]))
            self.fp_TOP.append(np.array([int(x) for x in rd_top.ToBitString()]))
            self.fp_MOR.append(np.array([int(x) for x in rd_mor.ToBitString()]))
#            print('folder {} was succeessfully operated !' .format(self.folder_name[b]))
            
    def get_MAC(self):
        return np.array(self.fp_MAC)
    
    def get_TOP(self):
        return np.array(self.fp_TOP)
    
    def get_MOR(self):
        return np.array(self.fp_MOR)
    


# In[3]:


if __name__ == '__main__':
    my_class = gen_fingerprint_array('/home/petitcloud/electrolyte_air_screen')
    fp = my_class.get_MOR()
    for c in range(len(fp)):
        print(len(fp[c]))

