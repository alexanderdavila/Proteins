from re import S, split
from Bio import SeqIO
from Bio.PDB import PDBList
import numpy as np
import time
from collections import defaultdict
from itertools import *

n=17
numAdd=(n // 2)

'''reading id's from pdb'''
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

'''Selecting structures from PDB'''
def download_pdb(pdblist):
    pdbl = PDBList()

    for i in pdblist:
         pdbl.retrieve_pdb_file(i,pdir='PDB',file_format="pdb") #-> str
         #print(size(pdbl))

'''read sequence from pdb'''
def read_sequence(pdblist):
    for list in pdblist:
        PDBFile = "./PDB/pdb"+list+".ent"
        with open(PDBFile, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                pdb_file.close()
                return record.seq
    

'''Crea matriz nxn de ceros'''
def crear_matriz():
    array=np.zeros((21,21))
    return array


'''Crear diccionario'''
def create_dictionary(aminoacids):
    dictionary=defaultdict(lambda: 20,zip(aminoacids,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]))
    return dictionary


'''recorrer-secuencia distancia n'''
def recorrer_sequence(sequence,init,n):
    var=sequence[init]
    
    tupla=[]
    for i in range(0,len(sequence)-n):
        if sequence[i]==var:
            tupla.append([sequence[i],sequence[i+n]])

    return tupla

'''obtener coordenadas'''
def get_coordenadas(tupla,dictionary):
    item=[]
    for i in range(0,len(tupla)):
        for j in range(0,len(tupla[0])):
            item.append(dictionary[tupla[i][j]])
    
    return item

def llenar_matriz(tupla,matriz):
    init1=0
    init2=1
    x=tupla[init1]
    y=tupla[init2]
    for i in range(0,int(len(tupla)/2)):
        if matriz[x][y]==0:
            matriz[x][y]=1
            init1+=2
            init2+=2
        else:
            matriz[x][y]+=1
            init1+=2
            init2+=2
        if init2 <=len(tupla):
            x=tupla[init1]
            y=tupla[init2]
    return matriz

def element_exist(list,index):

    if list.count(list[index])>=2:
        flag=True
    else:
        flag=False
    return flag

# def ventana(sequence,init,wsize):

#     end=init+wsize
#     substring=sequence[init:end]

#     return substring

def prepararSecuencia(sequence):
    complemento=list("-"*numAdd)
    sequence_complete=complemento+list(sequence)+complemento
    #print(aminoacidos)
    return sequence_complete

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    l=[]
    l=tee(iterable,n)
    for i in range(1,n):
      for a in range(i):
        next(l[i],None)
    return zip(*l)

def crear_matrices(ventana,diccionario):
            matriz1=crear_matriz()
            matriz2=crear_matriz()
            matriz3=crear_matriz()
            lista=[]
            matrices=[]
            #if i <= len(sequence_original)-n:

                #print("\nsecuencia: "+sequence)
            #f.write("\nsecuencia: "+str(ventana))
            for j in range(0,len(ventana)-1):
                lista.append(ventana[j])
                flag=element_exist(lista,j)
                if flag:
                    continue
                coordenadas1=recorrer_sequence(ventana,j,1)
                tupla1=get_coordenadas(coordenadas1,diccionario)
                matrizd1=llenar_matriz(tupla1,matriz1)
                #print(coordenadas1)
                #print("coordenadas1: "+ str(tupla1))
                if j < len(ventana)-2:
                    coordenadas2=recorrer_sequence(ventana,j,2)
                    tupla2=get_coordenadas(coordenadas2,diccionario)
                    matrizd2=llenar_matriz(tupla2,matriz2)
                    #print("coordenadas2: "+ str(tupla2))
                    #print(coordenadas2)
                if j < len(ventana)-3:
                    coordenadas3=recorrer_sequence(ventana,j,3)
                    tupla3=get_coordenadas(coordenadas3,diccionario)
                    matriz3d=llenar_matriz(tupla3,matriz3)
                    #print("coordenadas3: "+ str(tupla3))
                    #print(coordenadas3)
            matrices=[matriz1,matriz2,matriz3]
            return matrices

def get_iterator(sequence_original):
    ventaneo=pairwise(sequence_original)
    a=combinations(ventaneo,2)
    return a

def main():
     f = open('test.txt','w')
     inicio=time.time()
     pdblist=id_reading("pruebaId.txt")
     #download_pdb(pdblist)
     #wsize=17
     aminoacids="DERKNHQSTAGVPLFYIMWC-"
     diccionario=create_dictionary(aminoacids)
     #print(diccionario)
     f.write(str(diccionario))
     for list in pdblist:   
        PDBFile = "./PDB/pdb"+list+".ent"
        #print("./PDB/pdb"+list+".ent")
        f.write("./PDB/pdb"+list+".ent")
        with open(PDBFile, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                pdb_file.close()
        #print("secuencia: "+record.seq)
        f.write("secuencia: "+str(record.seq))
        #print("tamaño secuencia: "+str(len(record.seq)))
        f.write("tamaño secuencia: "+str(len(record.seq)))
        #sequence_original="--------"+read_sequence()+"--------"
        sequence_original=prepararSecuencia(record.seq)
        a=get_iterator(sequence_original)
        matricesvA=[]
        matricesvR=[]
        matrices={} 
        for vA,vR in a:
            matricesvA.append(crear_matrices(vA,diccionario)) 
            matricesvR.append(crear_matrices(vR,diccionario))
           
        

        fin=time.time()
        print("TIEMPO:"+str(fin-inicio))
     f.close()    
     
     matrices={"matricesvA":matricesvA,"matricesvR":matricesvR}
     print(len(matrices['matricesvR']))
if __name__ == '__main__':
    main()    


# '''Separa la secuencia'''
# def split_seq(sequence):
#     sequence=str(sequence)
#     #print(sequence)
    
#     split=[]
#     init=0
#     end=1
#     #print(len(sequence))
#     iterator=(len(sequence))
    
#     for i in range(int(iterator)):
#         split.append(sequence[init:end])
#         init=end
#         end=end+1    
#     return split


# def mezclar_lista(lista_original):
#     # Crear una copia, ya que no deberíamos modificar la original
#     lista = lista_original[:]
#     # Ciclo for desde 0 hasta la longitud de la lista -1
#     longitud_lista = len(lista)
#     for i in range(longitud_lista):
#         indice_aleatorio = random.randint(0, longitud_lista - 1)
#         # Intercambiar
#         temporal = lista[i]
#         lista[i] = lista[indice_aleatorio]
#         lista[indice_aleatorio] = temporal
#     # Regresarla
#     return lista

# def sonListasIguales(lista_a, lista_b):
 
#     if len(lista_a) != len(lista_b):
#         return False
#     else:
#         for i in range(0,len(lista_a)):
#             if(lista_a[i] != lista_b[i]):
#                 return False
 
#     return True
