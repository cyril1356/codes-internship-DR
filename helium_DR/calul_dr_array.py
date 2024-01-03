import csv
import math
from scipy import constants
import numpy as np
import time
 ###if we want to have a beautiful code we could work with object code.
t_0 = time.time()
def transform(
    mat, indices
): 
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
compt=0
while np.array_equal(np.where(a == Energie_d[compt][0])[0],np.array([])):## don't forget np.were give a tuple
    compt+=1
ind = np.where(a == Energie_d[compt][0])[0][
    0
]
##first indice of the d level (they are the bigger ones in DR.lev)
TR = TR[ind:]  ##we take only the d_level to be the upper states

TR = TR[TR[:, 2].argsort()]
a = TR.T[2]



## porblem we never have the last linal state in the lower level in the data base tr
compt=-1
while np.array_equal(np.where(a == Energie_f[compt][0])[0],np.array([])):## don't forget np.were give a tuple
    compt-=1
ind = np.where(a == Energie_f[compt][0])[0][
    0
]  ##take the last indice of the f level(they are the smaller ones on DR.lev)

TR = TR[: ind + 1]  ##we take only the f level to be the lower states
TR = TR[TR[:, 0].argsort()]

##np.savetxt('arrC.csv', TR[:,[0,2]], delimiter=',')

#last modification of tr we got a problem if there aren't any radiative decay d==>f were d is the fist(resp last) level indice

##ok sames indices
np.savetxt('TR_original+2p2d',TR[:,[0,2,-1]])

##we have the good tr now

##print(len(Energie_d),len(sum_same_indices(TR[:,[0,-1]])),len(sum_same_indices(AI[:,[0,-1]])))#ok =)
if TR[-1][0]!= Energie_d[-1][0]:
    a=np.zeros((1,len(TR[-1])))
    a[0][0]=Energie_d[-1][0]
    TR=np.concatenate((TR,a),axis=0)

if TR[0][0]!= Energie_d[0][0]:##problem here !!
    a=np.zeros((1,len(TR[-1])))
    a[0][0]=Energie_d[0][0]
    TR=np.concatenate((a,TR),axis=0)###not tcheked watch out 





Tr_sum = sum_same_indices(TR[:, [0, -1]])


##ok we got everything

##we only need of the value Aa[0] and Aa[2]
d_indices = Energie_d.T[0]


Ai = np.loadtxt(
    "Ai.csv", delimiter=","
)  # (resonnace energie, partial autoimization rate (d->i) and the total autonization rate)

d_indices = Energie_d.T[0]

'''Ai = np.zeros((len(d_indices), 3))
com = 0
for element in d_indices:
    Ai[com]=np.array(fac.AIBranch("DR.ai.b",int(element),int(Energie_i[0][0]))[:-1])
    com+=1'''

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
t_2 = time.time() - t_0
print(t_2, "temps calcul np")







np.savetxt("S_strontum.csv", S, delimiter=" ")

width = (Ai_2 + Ar_2) * 6.582119569e-16  ## reduced constant plank in ev

np.savetxt("S_strontum.csv", TR[:,[0,2,-1]], delimiter=" ")


results = np.stack((Ai_0, width, S), axis=-1)

ind=np.where(results[:,0]>10**(2))
print(ind)
#results=np.delete(results,ind,0)##

np.savetxt("E_W_S_strontum.csv", results, delimiter=" ")


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
        *np.power( math.e,-Ai_0 * constants.electron_volt / (T * constants.Boltzmann)))##boltazmann contant in J.K not eV.Kelvin 


    a = np.sum(a)
    factor = constants.h**3 / pow(
        (2 * constants.pi * constants.electron_mass * constants.Boltzmann * T), 3 / 2
    )
    #print(a,Ai_0,Ai_mat)
    return factor * a * (10**6)
    """      * (
            (math.e)
            ** (-inv_resonance * constants.electron_volt / (T * constants.Boltzmann))
        )"""


T = 1160450.500617  # for C
T = 11604500.500617  # for S

# [Aa[0],Bd,TwoJd,TwoJi,Aa[1],S])
matrix_ratecoef = [Ai_0, Bd, twojd, twoji, Ai_1, S]


a = np.linspace(20, 10**5, 500000)
rate = np.zeros(500000)
for i in range(len(rate)):
    rate[i] = calcul_DR_rates(matrix_ratecoef, a[i] * 11606)
print(len(rate))
# save arr to a CSV file using numpy.savetxt()
np.savetxt("DR_rates_strontum_2.csv", rate, delimiter=" ")
print("DR rate :", calcul_DR_rates(matrix_ratecoef, T))


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
    a = factor * np.power(math.e,a) * matrix_ratecoef[5]*constants.electron_volt#barn and ev
    return np.sum(a)


a = np.linspace((results[0][0]) *constants.electron_volt, (results[-1][0]*1.1)* constants.electron_volt, 50000)
rate = np.zeros((50000,2))
for i in range(len(rate)):
    rate[i][1] = cross_section(matrix_ratecoef, 2 * constants.electron_volt, a[i])
    rate[i][0]=a[i]/constants.electron_volt


np.savetxt("arr_cross_strontum_2.csv", rate, delimiter=" ")