"""
        EE2703 Assignment 1 submission by Aniket Kukreti EP19B036
"""

from sys import argv,exit

CIRCUIT = '.circuit'
END = '.end'

fname=argv[1]
if len(argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()
if '.netlist' not in fname:
    print('\nInvalid file extension. Please enter a .netlist file.')
    exit()

try:
    with open(fname) as f:
        lines = f.readlines()
        f.close()
        start = -1; end = -2
        for line in lines:              
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
                break
        if start >= end or start<0:                
            print('\nInvalid circuit definition.')
            exit()
        
        l2=lines[start+1:end]   #Checking branch definitions
        T1 = ["R","L","C","V","I"]
        T2 = ["H","F"]
        T3 = ["E","G"]
     
        for line in l2:
            if any(x in line.split('#')[0].split()[0] for x in T1):
                if len(line.split('#')[0].split()) != 4:
                    print('Invalid branch definition. Please check .netlist file.')
                    print("Line:"+line)
                    exit()
                    
            if any(x in line.split('#')[0].split()[0] for x in T2):
                if len(line.split('#')[0].split()) != 5:
                    print('Invalid branch definition. Please check .netlist file.')
                    print("Line:"+line)
                    exit()
                    
            if any(x in line.split('#')[0].split()[0] for x in T3):
                if len(line.split('#')[0].split()) != 6: 
                    print('Invalid branch definition. Please check .netlist file.')
                    print("Line:"+line)
                    exit()  
                
        for line in reversed([' '.join(reversed(line.split('#')[0].split())) for line in l2]):
            print(line)

except IOError:
    print('\nInvalid file. Please check the name of the input file and try again.')
    exit()



