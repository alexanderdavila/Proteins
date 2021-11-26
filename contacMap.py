import Bio.PDB
import numpy
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from numpy.core.fromnumeric import size
from numpy.lib.shape_base import array_split, split
import pylab
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import seaborn as sns

#from rpy2.robjects import r

# pdb_code = "4nkk"
# pdb_filename = "pdb4nkk.ent" #not the full cage

'''FUNCTIONS '''
def id_reading (file):
    fic = open(file, "r")
    lines = []
    for line in fic:
        lines.append(line)
    fic.close()
    PDBlist2=[]
    PDBlist2=lines[0].split(',')
    #print(PDBlist2)
    #print(size(PDBlist2))
    return PDBlist2

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    try:
        diff_vector = residue_one["CA"].coord -residue_two["CA"].coord
        distance=0
        distance=numpy.sqrt(numpy.sum(diff_vector*diff_vector))
        replace=numpy.nan_to_num(distance,nan=-1)
        #print(replace)
        return replace
    except:       
        return -1
   

def calc_dist_matrix(chain_one, chain_two) :
    #print (chain_one)
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

def create_contact_map(distance_matrix) :
    distance_matrix[distance_matrix>8]=0
    distance_matrix[distance_matrix>0]=1
    
    return distance_matrix

def pintar(contact_map,name):
    dct = {1: 5., 0: 1., -1: 10.}
    n = [[dct[i] for i in j] for j in contact_map]
    plt.imshow(n, cmap='brg', vmin=1, vmax=10)
    #print(contact_map)
    plt.savefig(name+'.png')
    #plt.show()
    

'''Selecting structures from PDB'''
def download_pdb(pdblist):
    pdbl = PDBList()

    for i in pdblist:
         pdbl.retrieve_pdb_file(i,pdir='PDB',file_format="pdb") #-> str
         #print(size(pdbl))
    
         

def create_dist_matrix(pdb_code):
    parser=PDBParser()
    structure = parser.get_structure(pdb_code, './PDB/pdb'+pdb_code+'.ent')
    model = structure[0]
    dist_matrix = calc_dist_matrix (model ["A"], model ["A"]) 
    return dist_matrix

def contact_map(pdblist):
    for i in pdblist:
        dist_matrix=create_dist_matrix(i)
        contact_map=create_contact_map(dist_matrix)
        pintar(contact_map,i)
        


def main():
    #guarda en una lista los Id's
    pdblist=id_reading("pruebaId.txt")
    #descarga los pdb's
    download_pdb(pdblist)
    #crea, pinta y guarda los mapas de contacto
    contact_map(pdblist)
    
    
main()



