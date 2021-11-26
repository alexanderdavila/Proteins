from itertools import *
from Bio import SeqIO
import numpy as np

'''read sequence from pdb'''
def read_sequence():
    PDBFile = "./PDB/pdb5jpj.ent"
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            #print('>' + record.id)
            #print(record.seq)
            pdb_file.close()
            return record.seq

def pairwise(iterable,n):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    l=[]
    l=tee(iterable,n)
    for i in range(1,n):
      for a in range(i):
        next(l[i],None)
    return list(zip(*l))

'''Separa la secuencia'''
def split_seq(sequence):
    sequence=str(sequence)
    print(sequence)
    
    split=[]
    init=0
    end=1
    #print(len(sequence))
    iterator=(len(sequence))
    
    for i in range(int(iterator)):
        split.append(sequence[init:end])
        init=end
        end=end+1    
    return split

'''Crea matriz nxn de ceros'''
def crear_matriz():
    array=np.zeros((21,21))
    return array


'''Crear diccionario'''
def create_dictionary(aminoacids):
    dictionary=dict(zip(aminoacids,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]))
    return dictionary

def recorrer_sequence(sequence,init,n):
    var=sequence[init]
    
    tupla=[]
    for i in range(0,len(sequence)-n):
        if sequence[i]==var:
            tupla.append([sequence[i],sequence[i+n]])

    return tupla

def main():
  aminoacids="DERKNHQSTAGVPLFYIMWC-"
  create_dictionary(aminoacids)
  sequence_completa=read_sequence()
  iterable=split_seq(sequence_completa)
  print(iterable)
  secuencia=pairwise(iterable,17)
  #print(secuencia[0])
  tupla=recorrer_sequence(secuencia[0],8,1)
  print("tupla:"+str(tupla))

if __name__ == '__main__':
    main()   