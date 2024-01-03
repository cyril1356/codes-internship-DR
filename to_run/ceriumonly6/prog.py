from pfac import fac
import csv
import math
from scipy import constants
import numpy as np
import time

# could create a file with everything
import sys
import os
import shutil

# IMPORTANT LOOK IF THE COMPTEUR IS AT 0 AT THE BEGUINNING OF THE CODE


##all functions


# Class of different styles
class style:
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    UNDERLINE = "\033[4m"
    RESET = "\033[0m"


def transform(mat, indices):
    """
    this function transform a matrice with string seraparte (givent by read_line appiled to the differnts files ) into
    a numpy matrice with float coefficent inside and return only the collums gave by the indices

    exemple
    transform([['2    2   3'],['3    3,4E-02     3']],[1 ,2])==> np.array([[2,3.40000e-2],[3,3]])

    """
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
    """
    this function work only for the files of types 1,2,4 of the data base of FAC
    this function take the file of one of the data base converted to a matrix by read lines
     and give only the lines of the results in the differents block as numpy matrices
    if we want to concatenate our differents blocks into 1 numpy matrix we have to say concatenate= True

    """
    nb_blocks = int(mat[6].split()[2])
    res = []
    data_base_type = int(mat[3].split()[2])
    T = False
    if data_base_type == 1:  ##energy
        start = 7  ## indice of the begining of the data base
        ##valuechange_number=3 no need always 3
        lenght_trans = 4  ## lenght of the transition between two blocks
        T = True  ##T is use to say that we won't have a grid
    if data_base_type == 2:  ## tr (radiative decays)
        start = 6
        lenght_trans = 6
        T = True
    if T:
        v = 0  ## will change the starting point
        ini = start + lenght_trans + 1
        for i in range(nb_blocks):
            numtrans = int(
                mat[start + v + 3].split()[2]
            )  ## number of element in the i^th block
            res.append(transform(mat[ini + v : ini + numtrans + v], indices))
            # mat[ini + v : ini + numtrans + v] is the matrice witch have the element of the i^th block we change it to a numpy matrice
            v += numtrans + lenght_trans

    if data_base_type == 5:
        start = 6
        lenght_trans = 5  ###watch out we have a grid
    if not T:
        v = 0
        ini = start + lenght_trans + 1
        for i in range(nb_blocks):
            grid = int(
                mat[start + v + 5].split()[2]
            )  ## nuber of element to caracterise the grid need we don't need this information
            numtrans = int(mat[start + v + 3].split()[2])
            res.append(
                transform(mat[ini + v + grid : ini + numtrans + v + grid], indices)
            )
            v += numtrans + lenght_trans + grid
    if concatenate:
        res = np.concatenate(tuple(res), axis=0)
    return res


def transform_range_i_f(
    mat, indices, ranges
):  ##supp indices i and f are on the 0,2  position
    ## not important in this code
    """
    transform a matrices but take only the incices witch are in the ranges values
    """
    range_i = ranges[0]
    range_f = ranges[1]
    r = np.zeros((len(mat), len(indices)))
    for i in range(len(mat)):
        res = np.zeros(len(indices))
        for j in range(len(indices)):
            val = mat[i].split()
            delel = []
            if (
                range_i[0] <= float(val[0]) <= range_i[1]
                and range_f[0] <= float(val[0]) <= range_f[1]
            ):
                delel.append(i)
            a = mat[i].split()[indices[j]]
            res[j] = float(a)
            r[i] = res

        r = np.delete(r, delel, axis=0)
    return r


def sum_same_indices(mat):
    ##have to be sorted by indices
    ##coulb be maybe done but creating sub matrices with the sames indices but I dont know how to do it with somthing OPTIMIZED!!
    ##this function will give the sum of the diffiertents indices
    ##this function can be use to extend our matrices too if they have different d_
    # costly in time use it after change TR into a smaller matrix with only the good rad decay
    # I will never outch it again

    """
    this function takes a matrices Nx2 soterd with respect to the first collumn and give a new numpy matrice
    were every lines who have the same fist value will be sum on the second colone

    exemple :
    sum_same_indices(np.array([[5,2],[7,3],[7,8],[12,8]]))
    return :
    np.array([[ 5.  2.]
    [ 6.  0.]
    [ 7. 11.]
    [ 8.  0.]
    [ 9.  0.]
    [10.  0.]
    [11.  0.]
    [12.  8.]])


    """
    d_ini = mat[0][0]
    d_end = mat[-1][0]
    num_d = d_end - d_ini + 1
    com1 = d_ini
    com2 = 0
    res = np.zeros((int(num_d), 2))
    sum = 0
    for i in range(len(mat)):
        if mat[i][0] == com1:
            sum += mat[i][1]
        else:
            res[com2] = np.array([com2 + d_ini, sum])
            com1 = mat[i][0]
            while com2 < com1 - d_ini:
                com2 += 1
                res[com2] = np.array([com2 + d_ini, 0])
            sum = mat[i][1]
    res[com2] = np.array([com2 + d_ini, sum])
    return res
    """
    could have use this !!
    more optimized !

    # create an numpy array
    a = np.array([1, 2, 3, 4, 8, 6, 2, 3, 9, 10])

    values = np.array([2, 3, 10])

    # index of first occurrence of each value
    sorter = np.argsort(a)

    print("index of first occurrence of each value: ",
        sorter[np.searchsorted(a, values, sorter=sorter)])
    """


##maybe usfull '1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6 6s2 4f14 5d10 6p6 7s2 5f14 6d10 7p6 8s2 5g18'


def give_electron_config_one_electron(closed_config, d_state):
    ##ok

    def concatenate(mat):
        res = ""
        for element in mat:
            element = element + " "
            res += element
        return res

    ###could creatie a progrma giving all the configuration callsed by n

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"
    config = config.split()
    closed_config = closed_config.split()

    indice_stop = config.index(d_state)
    indic_ini = config.index(closed_config[-1]) + 1
    res = []
    for i in range(indic_ini, indice_stop):
        res.append(concatenate(config[: i + 1]))

    return res


def give_electron_config(closed_config: str, d_state: str):
    """
    give all the 'f' states for a given ''d_state' where d_state is a str with the new shell of the bound
    electron and the new shell of the free electron
    exemple :
    give_electron_config('1s 2s 2p 3s 3p 4s 3d','5p 5d') ==>['1s 2s 2p 3s 3p 4s 3d 4p 5d', '1s 2s 2p 3s 3p 4s 3d 4p 5s 5d',
    '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5d', '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p ', '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s ',
    '1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f ']
    """

    def concatenate(mat):
        res = ""
        for element in mat:
            element = element + " "
            res += element
        return res

    def conca(mat):
        res = []
        for element in mat:
            res += element

    def ajoute_str(mat, state):
        ##change matrix watchout
        for i in range(len(mat)):
            mat[i] += state

        return mat

    res = []
    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"
    config = config.split()
    d_state = d_state.split()
    ##bound electron
    bound = d_state[0]
    indice_stop1 = config.index(
        bound
    )  ##subset include in each othe for the radiative decay

    res.append(
        ajoute_str(
            give_electron_config_one_electron(closed_config, config[indice_stop1]),
            d_state[1],
        )
    )

    ##free electron rad decay
    indice_stop2 = config.index(d_state[1])
    res.append(
        give_electron_config_one_electron(
            concatenate(config[:indice_stop1]), d_state[1]
        )
    )

    new_res = sum(res, [])
    return new_res


def give_electron_config_FAC(closed_config: str, d_state: str):
    """
    function made for fac calculations:
    wath out closed config no comma only spaces, closed config put the 'i' state inside too.

    exemple:
    give_electron_config_FAC('1s 2s 2p 3s 3p 4s 3d 4p 5s','5p 5g') ==> ['4d 5g', '5p 8s', '5p 7p', '5p 6d', '5p 5f', '5p 7s', '5p 6p', '5p 5d', '5p 4f', '5p 6s', '5p 5p']
    """

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"
    config = config.replace(closed_config, "")
    config = config.split()
    d_state = d_state.split()
    ##bound electron
    bound = d_state[0]
    indice_stop1 = config.index(bound)
    indice_ini = config.index(d_state[1])
    res = []
    for i in range(indice_stop1):
        a = config[i] + " " + d_state[1]
        res += [a]
    # free electron

    for i in range(
        indice_ini - 1, indice_stop1, -1
    ):  ##put -1 if you want have config 5p 5p for d_state[0]=5p
        res += [d_state[0] + " " + config[i]]
    return res


def give_electron_config_stars(closed_config, transition, end_stars):
    """
    gives all stars configurastion
    transition: explain the transition on witch the bound electron will go ex 5p
    list_states ground: are all the values configuration for n_ground: exemple for the transition 5p we have to put
    end star is an int saying the shell in witch the d electron will be catch
    """

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"
    config = config.split()
    closed_config = closed_config.split()
    indice_bound_transition = config.index(transition)
    value_n_bound = config[indice_bound_transition][0]
    list_bound = []

    # bound electrons
    for element in config:
        if element not in closed_config:
            if element[0] == value_n_bound:
                if element == config[indice_bound_transition]:
                    list_bound.append(transition + "2 ")
                else:
                    if config.index(element) > indice_bound_transition:
                        list_bound.append(
                            config[indice_bound_transition] + " " + element
                        )
                    else:
                        list_bound.append(
                            element + " " + config[indice_bound_transition]
                        )

    for i in range(int(value_n_bound) + 1, end_stars + 1):
        list_bound.append(str(i) + "* ")
    return list_bound





## need of a compteur so as to know in witch run we are

compteur = np.loadtxt("compteur.csv", delimiter=",")

compteur = int(compteur)


# start FAC calculations


t_0 = time.time()

# fac.InitializeMPI(64)


# It is important to include mixing of the initial ground state with
# excited states in the ground configuration (if any exist)
# we suppose taht there are no exitation and for Sr II !
# #1s2 	2s2 	2p6 	3s2 	3p6 	3d10 	4s2 	4p6 	4d10 	4f14 	5s2 	5p6 	5d10 	5f14 	6s2 	6p6 	6d10 	7s2 	7p6 8s2 5g18 6f14 7d10 8p6 9s2

# fac.Config('i', '1s2 2*8 3*18 4s2 4p6 5s')
# "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p     6s 4f 5d 6p 7s 5f 6d 7p 8s 5g 9s 8p 7d 6d"

# all_config_possible='1s 2s 2p 3s 3p 3d 4s 4p 4d 4f 5s 5p 5d 5f 5g 6s 6p 6d 6f 6g 6h 7s 7p 7d 7f 7g 7h 7i 8s 8p 8d 8f 8g 8h 8i 8k 9s 9p 9d 9f 9g 9h 9i 9k 9l34 10s2 10p6 10d10 10f14 10g18 10h22 10i26 10k30 10l34 10m38 11s2 11p6 11d10 11f14 11g18 11h22 11i26 11k30 11l34 11m38 11n42 12s2 12p6 12d10 12f14 12g18 12h22 12i26 12k30 12l34 12m38 12n42 12o46 13s2 13p6 13d10 13f14 13g18 13h22 13i26 13k30 13l34 13m38 13n42 13o46 13q50 14s2 14p6 14d10 14f14 14g18 14h22 14i26 14k30 14l34 14m38 14n42 14o46 14q50 14r54 15s2 15p6 15d10 15f14 15g18 15h22 15i26 15k30 15l34 15m38 15n42 15o46 15q50 15r54 15t58 16s2 16p6 16d10 16f14 16g18 16h22 16i26 16k30 16l34 16m38 16n42 16o46 16q50 16r54 16t58 16u62 17s2 17p6 17d10 17f14 17g18 17h22 17i26 17k30 17l34 17m38 17n42 17o46 17q50 17r54 17t58 17u62 17v66 18s2 18p6 18d10 18f14 18g18 18h22 18i26 18k30 18l34 18m38 18n42 18o46 18q50 18r54 18t58 18u62 18v66 18w70 19s2 19p6 19d10 19f14 19g18 19h22 19i26 19k30 19l34 19m38 19n42 19o46 19q50 19r54 19t58 19u62 19v66 19w70 19x74 20s2 20p6 20d10 20f14 20g18 20h22 20i26 20k30 20l34 20m38 20n42 20o46 20q50 20r54 20t58 20u62 20v66 20w70 20x74 20y78 21s2 21p6 21d10 21f14 21g18 21h22 21i26 21k30 21l34 21m38 21n42 21o46 21q50 21r54 21t58 21u62 21v66 21w70 21x74 21y78 21z82'


fac.SetAtom("Sr")
fac.Closed(
    " 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6"
)  # could the electron goses to the shell 4d ??## watch out here 4* directly

##config groung state cerium III    " 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 6s1"
fac.Config("i", "6s")

list_confi = [
    "6p 6p",
    "6p 7s",
    "6p 5f",
    "6p 6d",
    "6p 7p",
    "6p 8s",
    "6p 5g",
    "6p 9s",
    "6p 8p",
    "6p 7d",
    "6p 6d",
]


fac.Config("d", list_confi[compteur])

closed_config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p"

final_config = give_electron_config_FAC(closed_config, list_confi[compteur])

for element in final_config:
    fac.Config("f", element)
fac.Config("f", "6s2")


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

# we need to take the values of aionization and transform it the a matrix

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
    Ai[com] = np.array(fac.AIBranch("DR.ai.b", int(element), int(Energie_i[0][0]))[:-1])
    com += 1


np.savetxt("Ai_%s.csv" % (compteur), Ai, delimiter=",")
t_1 = time.time() - t_0


##part 2


filename = "DR.lev"
file0 = open(filename)

lines = file0.readlines()

a = transform_with_block(lines, [0, 2, 5])

# we take the indices , the energy in ev and the 2J values of the differtents configurations

Energie_i = a[1]
Energie_d = a[2]
Energie_f = a[0]

file0.close()

"""file0 = open("DR.ai")
lines = file0.readlines()


nb_block = int(lines[6].split()[2])

AI = transform_with_block(lines, [0, 1, 4, 5], True)##don't need it

AI = AI[AI[:, 0].argsort()]

file0.close()"""


Ai = np.loadtxt(
    "Ai_%s.csv" % (compteur), delimiter=","
)  # (resonnace energie, partial autoimization rate (d->i) and the total autonization rate)


file0 = open("DR.tr")

lines = file0.readlines()

TR = transform_with_block(lines, [0, 1, 2, 3, 4, 6], True)
##for radiative decay i==>f
# indices_i, 2J_i, indice_f,2J_f, E_i-E_f in ev, Tr rate 1/s,

TR = TR[TR[:, 0].argsort()]


# we sort our matrice with respect to the indices of the upper levels
file0.close()

##problem we can have for tr a radiative decay d==>d f==>f problem we only want d==>f
a = TR.T[0]
compt = 0
while np.array_equal(
    np.where(a == Energie_d[compt][0])[0], np.array([])
):  ## don't forget np.were give a tuple
    compt += 1
ind = np.where(a == Energie_d[compt][0])[0][0]
##first indice of the d level (they are the bigger ones in DR.lev)
TR = TR[ind:]  ##we take only the d_level to be the upper states

TR = TR[TR[:, 2].argsort()]
a = TR.T[2]

## porblem we never have the last linal state in the lower level in the data base tr
compt = -1
while np.array_equal(
    np.where(a == Energie_f[compt][0])[0], np.array([])
):  ## don't forget np.were give a tuple
    compt -= 1
ind = np.where(a == Energie_f[compt][0])[0][
    0
]  ##take the last indice of the f level(they are the smaller ones on DR.lev)

TR = TR[: ind + 1]  ##we take only the f level to be the lower states
TR = TR[TR[:, 0].argsort()]

##np.savetxt('arrC.csv', TR[:,[0,2]], delimiter=',')

# last modification of tr we got a problem if there aren't any radiative decay d==>f were d is the fist(resp last) level indice

##ok sames indices

##we have the good tr now

##print(len(Energie_d),len(sum_same_indices(TR[:,[0,-1]])),len(sum_same_indices(AI[:,[0,-1]])))#ok =)
if TR[-1][0] != Energie_d[-1][0]:
    a = np.zeros((1, len(TR[-1])))
    a[0][0] = Energie_d[-1][0]
    TR = np.concatenate((TR, a), axis=0)

if TR[0][0] != Energie_d[0][0]:  ##problem here !!
    a = np.zeros((1, len(TR[-1])))
    a[0][0] = Energie_d[0][0]
    TR = np.concatenate((a, TR), axis=0)  ###not tcheked watch out


Tr_sum = sum_same_indices(TR[:, [0, -1]])


##ok we got everything

##we only need of the value Aa[0] and Aa[2]
d_indices = Energie_d.T[0]


"""Ai = np.zeros((len(d_indices), 3))
com = 0
for element in d_indices:
    Ai[com]=np.array(fac.AIBranch("DR.ai.b",int(element),int(Energie_i[0][0]))[:-1])
    com+=1"""

##we got everything
assert (Energie_d[:, 0] == Tr_sum[:, 0]).all()
twojd = Energie_d[:, -1]
twoji = Energie_i[0][-1]
Ai_0 = Ai[:, 0]  # res energie in eV
Ai_1 = Ai[:, 1]  # partial autoinisation rate 1/s
Ai_2 = Ai[:, 2]  # total autoinization rate (d =>i)
Ar_2 = Tr_sum[:, 1]  # total radiative rate (d=>f)


inv_resonance_value = (
    2
    * pow(constants.hbar, 3)
    / (
        pow(
            Ai_0 * constants.electron_volt
            + constants.electron_mass * constants.c * constants.c,
            2,
        )
        / pow(constants.c, 2)  ##electron mass and photon velocity  and before h bar
        - pow(
            constants.electron_mass * constants.c, 2
        )  # dr resonance energy (look on internet) in energy unit  SOULD I MULT BY 2 ?? NOT THE RESONNANCE VALUE
    )
)

Bd = np.divide(
    Ar_2, Ai_2 + Ar_2, out=np.zeros_like(Ar_2), where=(Ai_2 + Ar_2) != 0
)  ## let say it is for the lauren gaussian
##return 0 if we want to divided by 0

S = (
    pow(constants.pi, 2)
    * inv_resonance_value
    * Ai_1  ## watch out it is not an enenergy
    * Bd
    * (twojd + 1)  # statistical weight
    / (twoji + 1)
    / 2
    / constants.electron_volt  ##return in ev
    * 10000  # cm^2
    / 1.0e-24  # to be in barn
)  # we make S finally

# Resonance strength in barn*eV:

# the result between my code and suvam code are a bit different because Ar_[2] (total radiative decay) calcul by fac is not the same (but almost the same) than the sum calculated here
##couls change that by creatind a matices like AI in the fac calcul

# Tabulate resonance energies, widths and strengths taking into accont only the autoinizations states possibles no need we got the sames results


np.savetxt("S_strontum_%s.csv" % (compteur), S, delimiter=" ")

width = (Ai_2 + Ar_2) * 6.582119569e-16  ## reduced constant plank in ev


results = np.stack((Ai_0, width, S), axis=-1)


ind = np.where(Ai[:, 0] > 10 ** (2))

if len(ind[0]) != 0:  ###false calculations
    print(
        ind,
        style.RED
        + "watch out there must have two sames d and f sates, or two sames f states in any case thoses values are not correct",
    )

    # create Directory (if needed)
    directory = "n=%s" % (compteur)

    # Parent Directory path
    # parent_dir = "/home/flury/Fileserver/test_big_calculs/"
    parent_dir = ""
    # Path
    path = os.path.join(parent_dir, directory)

    # Create the directory
    if not os.path.isdir("n=%s" % (compteur)):  ##seems to work but strange
        # os.mkdir(path)
        if not os.path.isdir(
            "n=%s__no_results_config__" % (compteur) + list_confi[compteur]
        ):
            os.mkdir("n=%s__no_results_config__" % (compteur) + list_confi[compteur])

    ##usless to put the fils in two if

    else:  # remove all data
        shutil.rmtree("n=%s" % (compteur))
        if not os.path.isdir(
            "n=%s__no_results_config__" % (compteur) + list_confi[compteur]
        ):
            os.mkdir("n=%s__no_results_config__" % (compteur) + list_confi[compteur])

    shutil.copyfile(
        parent_dir + "DR.lev",
        parent_dir
        + "n=%s__no_results_config__" % (compteur)
        + list_confi[compteur]
        + "/DR.lev",
    )  ##could put dr lev compteur too if we want
    os.remove(parent_dir + "DR.lev")

    shutil.copyfile(
        parent_dir + "DR.tr",
        parent_dir
        + "n=%s__no_results_config__" % (compteur)
        + list_confi[compteur]
        + "/DR.tr",
    )
    os.remove(parent_dir + "DR.tr")

    shutil.copyfile(
        parent_dir + "DR.ai",
        parent_dir
        + "n=%s__no_results_config__" % (compteur)
        + list_confi[compteur]
        + "/DR.ai",
    )
    os.remove(parent_dir + "DR.ai")

    score = give_electron_config_FAC(closed_config, list_confi[compteur])

    with open(
        "n=%s__no_results_config__" % (compteur) + list_confi[compteur] + "/config_f",
        "w",
    ) as f:
        for s in score:
            f.write(s + " ;  ")

    os.remove(parent_dir + "S_strontum_%s.csv" % (compteur))
    os.remove(parent_dir + "Ai_%s.csv" % (compteur))

    t_2 = time.time() - t_1 - t_0
    print("loop's number ", compteur, "finished")
    print("time FAC calculation :", t_1, "       time numpy calculs:", t_2)
    compteur += 1
    compteur = np.array([compteur])
    np.savetxt("compteur.csv", compteur, delimiter=",")

    print(
        style.GREEN
        + "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"
    )

    # reload the program
    if compteur < len(list_confi):
        os.execv(sys.executable, ["python"] + sys.argv)
    else:
        compteur = 0
        compteur = np.array([compteur])
        np.savetxt("compteur.csv", compteur, delimiter=",")
        quit()  ## end the programm

    # results=np.delete(results,ind,0)


np.savetxt("E_W_S_strontum_%s.csv" % (compteur), results, delimiter=" ")


def calcul_DR_rates(matrix_ratecoef, T):  ## works great
    """
    calcul of the resonnance strenth with the formula given in the dongs paper
    matrice_rate_coeff contains [Ai_0,Bd,TwoJd,TwoJi,Ai_1,S])

    """

    ##if it is too log can make into a matrix numpy
    ## remeber element =[inv_resonance_value,Bd,TwoJd,TwoJi,Aa[1]]
    gd = matrix_ratecoef[2] + 1
    gi = matrix_ratecoef[3] + 1  ###don't forget about the + 1
    Ai_mat = matrix_ratecoef[4]
    Ai_0 = matrix_ratecoef[0]
    Bd = matrix_ratecoef[1]
    a = (
        gd
        / (2 * gi)
        * Ai_mat
        * Bd
        * np.power(math.e, -Ai_0 * constants.electron_volt / (T * constants.Boltzmann))
    )  ##boltazmann contant in J.K not eV.Kelvin

    a = np.sum(a)
    factor = constants.h**3 / pow(
        (2 * constants.pi * constants.electron_mass * constants.Boltzmann * T), 3 / 2
    )
    # print(a,Ai_0,Ai_mat)
    return factor * a * (10**6)
    """      * (
            (math.e)
            ** (-inv_resonance * constants.electron_volt / (T * constants.Boltzmann))
        )"""


T = 1160450.500617  # for C
T = 11604500.500617  # for S

# [Aa[0],Bd,TwoJd,TwoJi,Aa[1],S])
matrix_ratecoef = [Ai_0, Bd, twojd, twoji, Ai_1, S]

taille = 5000
a = np.linspace(0.1, 5, taille)
b = np.linspace(5, 10**3, taille)
rate = np.zeros(taille * 2)
for i in range(len(rate)):
    if i < taille:
        rate[i] = calcul_DR_rates(matrix_ratecoef, a[i] * 11606)
    else:
        rate[i] = calcul_DR_rates(matrix_ratecoef, b[i - taille] * 11606)

# save arr to a CSV file using numpy.savetxt()
np.savetxt("DR_rates_strontum_%s.csv" % (compteur), rate, delimiter=" ")


## we will now calculate the DR cross section


def cross_section(matrix_ratecoef, delta_E, epsilon):
    """
    calul of the cross section by convolute eache resonnance strenghts with the resolution function (delta_E gaussian distribution)
    this formula is given in the dongs paper too
    """

    factor = 2 / delta_E * (math.log(2) / constants.pi) ** (1 / 2)
    a = (
        -4
        * math.log(2)
        * (((epsilon - matrix_ratecoef[0] * constants.electron_volt) / delta_E) ** 2)
    )
    a = (
        factor * np.power(math.e, a) * matrix_ratecoef[5] * constants.electron_volt
    )  # barn and ev
    return np.sum(a)


a = np.linspace(0.1 * constants.electron_volt, 5 * constants.electron_volt, 50000)
rate = np.zeros((50000, 2))
for i in range(len(rate)):
    rate[i][1] = cross_section(matrix_ratecoef, 0.05 * constants.electron_volt, a[i])
    rate[i][0] = a[i] / constants.electron_volt


np.savetxt("arr_cross_strontum_%s.csv" % (compteur), rate, delimiter=" ")

t_2 = time.time() - t_1 - t_0


# save data


# copy the contents of the demo.py file to  a new file called demo1.py


# use parent directory if you are working withought the cluster

# create Directory (if needed)
directory = "n=%s" % (compteur)

# Parent Directory path
# parent_dir = "/home/flury/Fileserver/test_big_calculs/"
parent_dir = ""
# Path
path = os.path.join(parent_dir, directory)

# Create the directory
if not os.path.isdir("n=%s" % (compteur)):  ##seems to work but strange
    # os.mkdir(path)
    os.mkdir("n=%s" % (compteur))


##copy usfull files

# parent_dir= '/home/flury/Fileserver/test_big_calculs/'
parent_dir = ""
shutil.copyfile(
    parent_dir + "Ai_%s.csv" % (compteur),
    parent_dir + "n=%s/Ai_%s.csv" % (compteur, compteur),
)
os.remove(parent_dir + "Ai_%s.csv" % (compteur))

shutil.copyfile(
    parent_dir + "arr_cross_strontum_%s.csv" % (compteur),
    parent_dir + "n=%s/arr_cross_strontum_%s.csv" % (compteur, compteur),
)
os.remove(parent_dir + "arr_cross_strontum_%s.csv" % (compteur))

shutil.copyfile(
    parent_dir + "DR_rates_strontum_%s.csv" % (compteur),
    parent_dir + "n=%s/DR_rates_strontum_%s.csv" % (compteur, compteur),
)
os.remove(parent_dir + "DR_rates_strontum_%s.csv" % (compteur))

shutil.copyfile(
    parent_dir + "E_W_S_strontum_%s.csv" % (compteur),
    parent_dir + "n=%s/E_W_S_strontum_%s.csv" % (compteur, compteur),
)
os.remove(parent_dir + "E_W_S_strontum_%s.csv" % (compteur))

shutil.copyfile(
    parent_dir + "DR.lev", parent_dir + "n=%s/DR.lev" % (compteur)
)  ##could put dr lev compteur too if we want
os.remove(parent_dir + "DR.lev")

shutil.copyfile(parent_dir + "DR.tr", parent_dir + "n=%s/DR.tr" % (compteur))
os.remove(parent_dir + "DR.tr")

shutil.copyfile(parent_dir + "DR.ai", parent_dir + "n=%s/DR.ai" % (compteur))
os.remove(parent_dir + "DR.ai")

os.remove(parent_dir + "S_strontum_%s.csv" % (compteur))
# could put somthing to delete the files

print("loop's number ", compteur, style.RED + "finished")
print(
    style.MAGENTA + "time FAC calculation :",
    t_1,
    style.MAGENTA + "      time numpy calculs:",
    t_2,
)
compteur += 1
compteur = np.array([compteur])
np.savetxt("compteur.csv", compteur, delimiter=",")

print(
    style.GREEN
    + "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"
)
print(style.WHITE + "   ")
# reload the program
if compteur < len(list_confi):
    os.execv(sys.executable, ["python"] + sys.argv)
else:
    compteur = 0
    compteur = np.array([compteur])
    np.savetxt("compteur.csv", compteur, delimiter=",")
