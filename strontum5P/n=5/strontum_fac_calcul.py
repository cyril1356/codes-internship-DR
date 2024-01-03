from pfac import fac
import csv
import math
from scipy import constants
import numpy as np 
import time
t_0=time.time()

fac.InitializeMPI(64)

def transform(
    mat, indices
):  ###if we want to have a beautiful code we could work with object code.
    r = np.zeros((len(mat), len(indices)))
    for i in range(len(mat)):
        res = np.zeros(len(indices))
        for j in range(len(indices)):
            a = mat[i].split()[indices[j]]
            res[j] = float(a)
        r[i] = res
    return r


def transform_with_block(
    mat, indices, concatenate=False
):  ##works only for energy , tr, ai
    nb_blocks = int(mat[6].split()[2])
    res = []
    data_base_type = int(mat[3].split()[2])
    T = False
    if data_base_type == 1:  ##energy
        start = 7
        ##valuechange_number=3 no need alway 3
        lenght_trans = 4
        T = True
    if data_base_type == 2:
        start = 6
        lenght_trans = 6
        T = True
    if T:
        v = 0
        ini = start + lenght_trans + 1
        for i in range(nb_blocks):
            numtrans = int(mat[start + v + 3].split()[2])
            res.append(transform(mat[ini + v : ini + numtrans + v], indices))
            v += numtrans + lenght_trans

    if data_base_type == 5:
        start = 6
        lenght_trans = 5  ###watch out we have a grid
    if not T:
        v = 0
        ini = start + lenght_trans + 1
        for i in range(nb_blocks):
            grid = int(mat[start + v + 5].split()[2])
            numtrans = int(mat[start + v + 3].split()[2])
            res.append(
                transform(mat[ini + v + grid : ini + numtrans + v + grid], indices)
            )
            v += numtrans + lenght_trans + grid
    if concatenate:
        res = np.concatenate(tuple(res), axis=0)
    return res



def create_config(closed_confi,d_state):


    def concatenate(mat):
        res=''
        for element in mat : 
            element=element + ' '
            res+=element
        return res

    config='1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g'
    config=config.split()
    closed_confi=closed_confi.split()
    indice_stop=config.index(d_state)
    indic_ini=config.index(closed_confi[-1])+1
    res=[]
    for i in range (indic_ini,indice_stop):
        res.append(concatenate(config[:i]))

    return res

fac.SetAtom('Sr')
fac.Closed('1* 2* 3* 4s2 4p6') #could the electron goses to the shell 4d ??## watch out here 4* directly



# It is important to include mixing of the initial ground state with
# excited states in the ground configuration (if any exist)
# we suppose taht there are no exitation and for Sr II !
# 1s, 2s, 2p 3s, 3p,4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p energy oes up

#fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')

#fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')
fac.Config('i', '5s')

fac.Config('d','5p2')
fac.Config('d','5p 4f')
fac.Config('d','5p 5d')
fac.Config('d','5p 5f')
fac.Config('d','5p 5g')




fac.Config('f','5s2')

fac.Config('f','5s 5p')
fac.Config('f','5s 4f')
fac.Config('f','5s 5d')
fac.Config('f','5s 5f')
fac.Config('f','5s 5g')


fac.Config('f','4d 5p')
fac.Config('f','4d 4f')
fac.Config('f','4d 5d')
fac.Config('f','4d 5f')
fac.Config('f','4d 5g')






##more exaustive list f

 

# Run FAC calculations for initial (ground) state energies

##calcul of the energy of the different configuration and construction of the potential
fac.SetBreit(-1, 2)
fac.ConfigEnergy(0)
fac.OptimizeRadial(["i"])
fac.ConfigEnergy(1)

# diagolalization of the hamiltonian and save the energyt in the most important data base
fac.Structure("DR.lev.b", ["i", "d", "f"])
fac.MemENTable("DR.lev.b")
fac.PrintTable("DR.lev.b", "DR.lev", 1)

# Calculate the Auger (autoionization) rates
fac.AITable(
    "DR.ai.b", ["d"], ["i"]
)  ##why do we have differents autoionization rates like our electron will always goes to d==>i (maybe different configuration inside the group d)
fac.PrintTable("DR.ai.b", "DR.ai", 1)

# Calculate the radiative decay rates


fac.TransitionTable("DR.tr.b", ["d"], ["f"])


fac.PrintTable("DR.tr.b", "DR.tr", 1)

#we need to take the values of aionization and transform it the a matrix 

filename = "DR.lev"
file0 = open(filename)

lines = file0.readlines()

a = transform_with_block(lines, [0, 2, 5])

Energie_i = a[1]
Energie_d = a[2]
Energie_f = a[0]

file0.close()

d_indices = Energie_d.T[0]

Ai = np.zeros((len(d_indices), 3))
com = 0
for element in d_indices:
    Ai[com]=np.array(fac.AIBranch("DR.ai.b",int(element),int(Energie_i[0][0]))[:-1])
    com+=1


np.savetxt('Ai.csv', Ai, delimiter=',')
t_1 = time.time() - t_0
print(t_1, "temps 1")



fac.ClearLevelTable()
fac.ClearOrbitalTable(0)



fac.SetAtom('Sr')
fac.Closed('1* 2* 3* 4s2 4p6') #could the electron goses to the shell 4d ??## watch out here 4* directly



# It is important to include mixing of the initial ground state with
# excited states in the ground configuration (if any exist)
# we suppose taht there are no exitation and for Sr II !
# 1s, 2s, 2p 3s, 3p,4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p energy oes up

#fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')

#fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')
fac.Config('i', '5s')

catch_shell=7
fac.Config('d','5p1 %s*1'%(catch_shell))

#rad decay of the free electron 
fac.Config('f','5p1 5g1')
fac.Config('f','5p1 5f1')
fac.Config('f','5p1 5d1')
fac.Config('f','5p1 4f1')

for i in range (6,catch_shell):
    fac.Config('f','5p1 %s*1'%(i))

#rad decay of the bound electron 
fac.Config('f','5s1 %s*1'%(catch_shell))
fac.Config('f','4d1 %s*1'%(catch_shell))





##more exaustive list f

 

# Run FAC calculations for initial (ground) state energies

##calcul of the energy of the different configuration and construction of the potential
fac.SetBreit(-1, 2)
fac.ConfigEnergy(0)
fac.OptimizeRadial(["i"])
fac.ConfigEnergy(1)

# diagolalization of the hamiltonian and save the energyt in the most important data base
fac.Structure("DR.lev.b", ["i", "d", "f"])
fac.MemENTable("DR.lev.b")
fac.PrintTable("DR.lev.b", "DR.lev", 1)

# Calculate the Auger (autoionization) rates
fac.AITable(
    "DR.ai.b", ["d"], ["i"]
)  ##why do we have differents autoionization rates like our electron will always goes to d==>i (maybe different configuration inside the group d)
fac.PrintTable("DR.ai.b", "DR.ai", 1)

# Calculate the radiative decay rates


fac.TransitionTable("DR.tr.b", ["d"], ["f"])


fac.PrintTable("DR.tr.b", "DR.tr", 1)

#we need to take the values of aionization and transform it the a matrix 

filename = "DR.lev"
file0 = open(filename)

lines = file0.readlines()

a = transform_with_block(lines, [0, 2, 5])

Energie_i = a[1]
Energie_d = a[2]
Energie_f = a[0]

file0.close()

d_indices = Energie_d.T[0]

Ai = np.zeros((len(d_indices), 3))
com = 0
for element in d_indices:
    Ai[com]=np.array(fac.AIBranch("DR.ai.b",int(element),int(Energie_i[0][0]))[:-1])
    com+=1


np.savetxt('Ai.csv', Ai, delimiter=',')
t_1 = time.time() - t_0
print(t_1, "temps 1")




fac.FinalizeMPI()



'''


fac.Closed('1* 2* 3* 4s2 4p6') #could the electron goses to the shell 4d ??## watch out here 4* directly



# It is important to include mixing of the initial ground state with
# excited states in the ground configuration (if any exist)
# we suppose taht there are no exitation and for Sr II !


#fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')
fac.Config('i', '5s')


# Intermediate- and final-state configs
# we got tree electron in the O shell and in this config the electron of the O shell don't move


fac.Config('d', '5p 5p')
fac.Config('d', '5p 5d')
fac.Config('d', '5p 5f')
fac.Config('d', '5p 5g')
fac.Config('d', '5p 6*')
fac.Config('d', '5p 7*')

fac.Config('d', '5d 5d')
fac.Config('d', '5d 5f')
fac.Config('d', '5d 5g')
fac.Config('d', '5d 6*')
fac.Config('d', '5d 7*')


fac.Config('d', '5f 5f')
fac.Config('d', '5f 5g')
fac.Config('d', '5f 6*')
fac.Config('d', '5f 7*')



fac.Config('d', '5g 5g')
fac.Config('d', '5g 6*')
fac.Config('d', '5g 7*')



fac.Config('f', '5* 5p')
fac.Config('f', '5* 5d')
fac.Config('f', '5* 5f')
fac.Config('f', '5* 5g')
fac.Config('f', '5* 6*') # really rare
fac.Config('f', '5* 7*')



'''



'''for i in range (6,8):
    fac.Config('d','5d %s*1'%(i))# the % signe say that we change the %s by i
    fac.Config('d','5p %s*1'%(i))
    fac.Config('d','5f %s*1'%(i))
    fac.Config('d','5g %s*1'%(i))


for i in range (6,8):    
    for j in range (i,8):
        fac.Config('d','%s*1 %r*1'%(i,j+1))
'''




"""
d 1s2 2*8 3*18 4s2 4p6  6*1 7*1
d 1s2 2*8 3*18 4s2 4p6  6*1 8*1
d 1s2 2*8 3*18 4s2 4p6  6*1 9*1
d 1s2 2*8 3*18 4s2 4p6  6*1 6*1
d 1s2 2*8 3*18 4s2 4p6  7*1 8*1
d 1s2 2*8 3*18 4s2 4p6  7*1 9*1
d 1s2 2*8 3*18 4s2 4p6  7*1 7*1
d 1s2 2*8 3*18 4s2 4p6  8*1 9*1
d 1s2 2*8 3*18 4s2 4p6  8*1 8*1
"""

'''for i in range (6-1,8-1):    
    for j in range (i,8-1):
        fac.Config('f','%s*1 %r*1'%(i,j+1))
    fac.Config('f',' %s*1 %r*1'%(i,i))'''



'''

different configuration but the program dosen't work if i right like this 
f 1s2 2*8 3*18 4s2 4p6  5*1 6*1
f 1s2 2*8 3*18 4s2 4p6  5*1 7*1
f 1s2 2*8 3*18 4s2 4p6  5*1 8*1
f 1s2 2*8 3*18 4s2 4p6  5*1 5*1
f 1s2 2*8 3*18 4s2 4p6  6*1 7*1
f 1s2 2*8 3*18 4s2 4p6  6*1 8*1
f 1s2 2*8 3*18 4s2 4p6  6*1 6*1
f 1s2 2*8 3*18 4s2 4p6  7*1 8*1
f 1s2 2*8 3*18 4s2 4p6  7*1 7*1

'''

# Run FAC calculations for initial (ground) state energies

