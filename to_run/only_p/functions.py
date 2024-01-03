import os

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

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g"
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
    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g"
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

    '''
    function made for fac calculations: 
    wath out closed config no comma only spaces, closed config put the 'i' state inside too.

    exemple:
    give_electron_config_FAC('1s 2s 2p 3s 3p 4s 3d 4p 5s','5p 5g') ==> ['4d 5g', '5p 8s', '5p 7p', '5p 6d', '5p 5f', '5p 7s', '5p 6p', '5p 5d', '5p 4f', '5p 6s', '5p 5p']
    '''

    config = "1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p 8s 5g"
    config=config.replace(closed_config,'')
    config=config.split()
    d_state=d_state.split()
    ##bound electron
    bound = d_state[0]
    indice_stop1 = config.index(bound)
    indice_ini= config.index(d_state[1])
    res=[]
    for i in range (indice_stop1):
        a=(config[i]+' '+d_state[1])
        res+=[a]
    #free electron

    for i in range (indice_ini-1,indice_stop1,-1):
        res+=[d_state[0]+' '+config[i]]
    return res










'''
# Directory 
directory = "n=1"
  
# Parent Directory path 
parent_dir = "/home/flury/Fileserver/test_big_calculs/"
  
# Path 
path = os.path.join(parent_dir, directory) 
  
# Create the directory 
# 'GeeksForGeeks' in 
# '/home / User / Documents' 
os.mkdir(path) 




import shutil

# copy the contents of the demo.py file to  a new file called demo1.py
shutil.copyfile('DR.ai', 'n=0/jbbhud')
compteur=0

'''

# 1s, 2s, 2p 3s, 3p,4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p energy oes up



score=give_electron_config_FAC('1s 2s 2p 3s 3p 4s 3d 4p 5s','5p 4f')

# Class of different styles
class style():
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'

print(style.YELLOW + "Hello, World!")