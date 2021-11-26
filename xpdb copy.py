from itertools import  combinations_with_replacement
cont=0
file = open("combinaciones.txt","w")
for i in combinations_with_replacement('DERKNHQSTAGVPLFYIMWC',17): 
  print(i, end=' ')
  file.write(str(i))
  cont +=1 # aa ab ac ad bb bc bd cc cd dd

print("TOTAL COMBINACIONES:", cont)
file.close()