"""
        EE2703 Assignment 2 by Aniket Kukreti EP19B036
"""

from sys import argv,exit
import math
import numpy as np
import cmath as cm

CIRCUIT = '.circuit'
END = '.end'
GND= 'GND'
ac=0

#Defining component classes
class R:
    def __init__(self,name,n1,n2,value):
        self.name=name
        self.n1=n1
        self.n2=n2
        self.value=value
        def __repr__(self):
            return(self.name)
class L:
    def __init__(self,name,n1,n2,value):
        self.name=name
        self.n1=n1
        self.n2=n2
        self.value=value
        def __repr__(self):
            return(self.name)
class C:
    def __init__(self,name,n1,n2,value):
        self.name=name
        self.n1=n1
        self.n2=n2
        self.value=value
        def __repr__(self):
            return(self.name)
class V:
    def __init__(self,name,n1,n2,typ,value):
        self.name=name
        self.n1=n1
        self.n2=n2
        self.typ=typ
        self.value=value
        def __repr__(self):
            return(self.name)
class I:
    def __init__(self,name,n1,n2,typ,value):
        self.name=name
        self.n1=n1
        self.n2=n2
        self.typ=typ
        self.value=value
        def __repr__(self):
            return(self.name)
        
#Opening the netlist file               
fname=argv[1]
if len(argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()
if '.netlist' not in fname:
    print('\nInvalid file extension. Please enter a .netlist file.')
    exit()
    
#Reading the netlist file
try:
    with open(fname) as f:
        lines = f.readlines()
        f.close()
        start = -1; end = -2
        
        for line in lines:              
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                if end>0:
                    continue
                end = lines.index(line)
            elif '.ac' == line[:3]:
                ac=1
                freq=float(line.split('#')[0].split()[2])
                break
        
        if start >= end or start<0:                
            print('\nInvalid circuit definition.')
            exit()
except IOError:
    print('\nInvalid file. Please check the name of the input file and try again.')
    exit()        

#Storing the lines into lists     
l2=lines[start+1:end]
tkn=[]
for line in l2:
    tkn.append(line.split('#')[0].split())
comps=[[] for x in range(5)] #[[R],[L],[C],[V],[I]]
nodes=[]

#Storing the components and their token values.
try:
    if ac == 0:
        for i in range(len(l2)):
            nodes.append(tkn[i][1])
            nodes.append(tkn[i][2])
            if 'R' in tkn[i][0]:
                comps[0].append(R(tkn[i][0],tkn[i][1],tkn[i][2],(float)(tkn[i][3])))
            elif 'L' in tkn[i][0]:
                comps[1].append(L(tkn[i][0],tkn[i][1],tkn[i][2],(float)(tkn[i][3])))
            elif 'C' in tkn[i][0]:
                comps[2].append(C(tkn[i][0],tkn[i][1],tkn[i][2],(float)(tkn[i][3])))
            elif 'V' in tkn[i][0]:
                if tkn[i][3] == "ac":
                    print('Error in netlist. AC source detected but not declared.')
                comps[3].append(V(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],float(tkn[i][4])))
            elif 'I' in tkn[i][0]:
                comps[4].append(I(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],float(tkn[i][4])))
    else:
        for i in range(len(l2)):
            nodes.append(tkn[i][1])
            nodes.append(tkn[i][2])
            if 'R' in tkn[i][0]:
                comps[0].append(R(tkn[i][0],tkn[i][1],tkn[i][2],complex(tkn[i][3])))
            elif 'L' in tkn[i][0]:
                comps[1].append(L(tkn[i][0],tkn[i][1],tkn[i][2],float(tkn[i][3])))
            elif 'C' in tkn[i][0]:
                comps[2].append(C(tkn[i][0],tkn[i][1],tkn[i][2],float(tkn[i][3])))
            elif 'V' in tkn[i][0]:
                if tkn[i][3] == "dc":
                    comps[3].append(V(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],(float)(tkn[i][4])))
                else:
                    ph=float(tkn[i][5])
                    val=complex(float(tkn[i][4])/(2*math.sqrt(2)))*complex(math.cos(ph*math.pi/180),math.sin(ph*math.pi/180))
                    comps[3].append(V(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],val))
            elif 'I' in tkn[i][0]:
                if tkn[i][3] == "dc":
                    comps[4].append(I(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],float(tkn[i][4])))
                else:
                    ph=float(tkn[i][5])
                    val=complex(float(tkn[i][4])/(2*math.sqrt(2)))*complex(math.cos(ph*math.pi/180),math.sin(ph*math.pi/180))
                    comps[3].append(I(tkn[i][0],tkn[i][1],tkn[i][2],tkn[i][3],val))
except IndexError:
    print("Incorrect component definition. Please check netlist file. (indices)")
except ValueError:
    print("Incorrect component value. Please check netlist file. (values)") 

#Creating list of distinct nodes
nodes1=list(dict.fromkeys(nodes));
if GND not in nodes1:
    print('GND node not present. Please check the netlist.')
nodes1.sort();
nodes1.remove(GND)
nodes1.insert(0,GND)
nodeindex=dict.fromkeys(nodes1); 
i=0
for x in nodes1: #Associating dictionary of nodes with integer values
    nodeindex[x]=i
    i+=1

#Creating the structure of conductance matrix and b matrix
n=len(nodeindex); 
k=len(comps[3]); 
matM=np.zeros((n+k,n+k), np.complex)
matB=np.zeros(n+k, np.complex)
matM[0][0]=1 #Since GND node is always V0=0
#Contribution of resistors towards matrix M
for r in comps[0]:
    if r.n1 != GND:
        matM[nodeindex[r.n1]][nodeindex[r.n1]]+= 1/r.value
        matM[nodeindex[r.n1]][nodeindex[r.n2]]-= 1/r.value
    if r.n2 != GND:
        matM[nodeindex[r.n2]][nodeindex[r.n2]]+= 1/r.value
        matM[nodeindex[r.n2]][nodeindex[r.n1]]-= 1/r.value
#Contribution of inductors towards matrix M
for l in comps[1]:
    if ac==0:
        if l.n1 != GND: 
            matM[nodeindex[l.n1]][nodeindex[l.n1]]+= 1e250 #Inductor has infinite conductance in a DC circuit. 
            matM[nodeindex[l.n1]][nodeindex[l.n2]]-= 1e250 #Hence, large value of conductance is used to mimic infinity
        if l.n2 != GND:
           matM[nodeindex[l.n2]][nodeindex[l.n2]]+= 1e250
           matM[nodeindex[l.n2]][nodeindex[l.n1]]-= 1e250
    else:
        if l.n1 != GND:
            matM[nodeindex[l.n1]][nodeindex[l.n1]]+= 1/complex(0,2*math.pi*freq*l.value)
            matM[nodeindex[l.n1]][nodeindex[l.n2]]-= 1/complex(0,2*math.pi*freq*l.value)
        if l.n2 != GND:
           matM[nodeindex[l.n2]][nodeindex[l.n2]]+= 1/complex(0,2*math.pi*freq*l.value)
           matM[nodeindex[l.n2]][nodeindex[l.n1]]-= 1/complex(0,2*math.pi*freq*l.value)
#Contribution of capacitors towards matrix M
for c in comps[2]:
    if ac==1: #Note that for DC, capacitor is open circuited, so it doesn't contribute to KCL at that node.       
        if c.n1 != GND:
            matM[nodeindex[c.n1]][nodeindex[c.n1]]+=1/complex(0,-1/(2*math.pi*freq*c.value))
            matM[nodeindex[c.n1]][nodeindex[c.n2]]-=1/complex(0,-1/(2*math.pi*freq*c.value))
        if c.n2 != GND:
           matM[nodeindex[c.n2]][nodeindex[c.n2]]+=1/complex(0,-1/(2*math.pi*freq*c.value))
           matM[nodeindex[c.n2]][nodeindex[c.n1]]-=1/complex(0,-1/(2*math.pi*freq*c.value))
#Contribution of voltage sources towards matrix M and b
s=0
for v in comps[3]:
    matM[n+s][nodeindex[v.n1]]+=1
    matM[n+s][nodeindex[v.n2]]-=1
    if v.n1 != GND:
        matM[nodeindex[v.n1]][n+s]+=1
    if v.n2 != GND:
        matM[nodeindex[v.n2]][n+s]-=1
    matB[n+s]+=v.value
    s+=1
#Contribution of current sources towards matrix b
for i in comps[4]:
    if i.n1 != GND:
        matB[nodeindex[i.n1]]-=i.value
    if i.n2 != GND:
        matB[nodeindex[i.n2]]+=i.value

#Solving M and b
try:
    matX=np.linalg.solve(matM,matB)
except np.linalg.LinAlgError:
    print('Error while solving the circuit. Please check if circuit has been defined correctly in the .netlist file.')

#Correcting negligible values to zeros
for i in range(n+k):
    if abs(matX[i]) < 1e-12:
        matX[i]=complex(0,0)

#Printing node voltages and voltage source currents
if ac == 0:    
    for y in range(n):
        print('Voltage at node '+str(y)+': '+str(matX[y].real)+' V')
    for y in range(n,n+k):
        print('Current through '+comps[3][y-n].name+': '+str(matX[y].real)+' A')
else:
    for y in range(n):
        print('Voltage at node '+str(y)+': '+str(cm.polar(matX[y]))+' V_rms')
    for y in range(n,n+k):
        print('Current through '+comps[3][y-n].name+': '+str(cm.polar(matX[y]))+' A_rms')