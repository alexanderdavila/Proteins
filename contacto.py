import Bio.PDB
from Bio import SeqIO
import numpy
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from numpy.core.fromnumeric import size
from numpy.lib.shape_base import array_split, split
import pylab
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
#import seaborn as sns
from itertools import *
import time

inicio = time.time()
time.sleep(1)
def pairwise(iterable,n):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    l=[]
    l=tee(iterable,n)
    for i in range(1,n):
      for a in range(i):
        next(l[i],None)
    return zip(*l)

test_list = ['-','-','-','-','-','-','-','-','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','-','-','-','-','-','-','-','-']
l=pairwise(test_list,17)
#print(l)
fin = time.time()
#print(fin-inicio)

n=17
numAdd=(n // 2)
'''FUNCTIONS '''
def id_reading (file):
    fic = open(file, "r")
    lines = []
    for line in fic:
        lines.append(line)
    fic.close()
    PDBlist2=[]
    PDBlist2=lines[0].split(',')
    return PDBlist2

'''read sequence from pdb'''
def read_sequence():
    PDBFile = "./PDB/pdb7A58.ent"
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            return record.seq

def prepararSecuencia(aminoacido):
    complemento=list("-"*numAdd)
    aminoacidos=complemento+list(aminoacido)+complemento
    #print(aminoacidos)
    return aminoacidos

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    l=[]
    l=tee(iterable,n)
    for i in range(1,n):
      for j in range(i):
        next(l[i],None)
    return zip(*l)

def contact_map(pdblist):
    for i in pdblist:
        #print("MATRIZ DE DISTANCIA "+i)
        dist_matrix=create_dist_matrix(i)
        #print(dist_matrix)
        contact_map=create_contact_map(dist_matrix)
        
        return contact_map
        #pintar(contact_map,i)

def create_dist_matrix(pdb_code):
    parser=PDBParser()
    structure = parser.get_structure(pdb_code, './PDB/pdb7A58.ent')
    model = structure[0]
    dist_matrix = calc_dist_matrix (model ["A"], model ["A"]) 
    return dist_matrix

def calc_dist_matrix(chain_one, chain_two) :
    #print (chain_one)
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

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

def create_contact_map(distance_matrix) :
    distance_matrix[distance_matrix>8]=0
    distance_matrix[distance_matrix>0]=1
    return distance_matrix
    
def hayContacto(ventanaA,ventanaR,contacMap):
    l=[]
    inicio=0
    for i in range(len(contacMap)):
        next(ventanaA)
        for j in range(len(contacMap)):
            next(ventanaR)
            if(contacMap[i][j]==1):
                l=list(next(ventanaA))
                print(l)
                print(list(next(ventanaR)))
                print(contacMap[i][j])                
    return l

def main():
    #Obtener secuencia
    secuence=read_sequence()
    
    #Preparacion de secuencia, se agregan "-" en cada extremo de la sec
    #para poder generar correctamente las ventanas n aminoacidos    
    aminoPreparados=prepararSecuencia(secuence)
    ventaneo=pairwise(aminoPreparados)
    a,b=tee(ventaneo,2)
    #print(next(a))
    #print(next(b))
    #print(next(b))
    #Obtener Matriz de distancia
    pdblist=id_reading("pruebaID.txt")
    
    #generarMatrizDeContacto(aminoPreparados)  
    #guarda en una lista los Id's
    #Recorrer lista de ventanas y generar matriz de contacto
    #print(pdblist)
    c=contact_map(pdblist)
    hayContacto(a,b,c)
main()