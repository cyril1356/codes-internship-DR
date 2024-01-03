import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


print('ede')
list_confi=['5p 5p','5p 6s','5p 4f','5p 5d','5p 6p','5p 7s','5p 5f','5p 6d','5p 7p', '5p 8s', '5p 5g']


convolution=[]
taille=5000

a = np.linspace(0.1, 5, taille)
b=np.linspace(5,10**3,taille)
x=np.concatenate((a,b))

for i in range (1,len(list_confi)):

    y_1 = np.loadtxt('/home/flury/Fileserver/test_big_calculs_convolution/n=%s/DR_rates_strontum_%s.csv'%(i,i), delimiter=' ')
    

    plt.plot(x, y_1,label="configuration %s"%(list_confi[i]))



plt.legend(fontsize=20)

plt.yscale('log')
plt.xscale('log')
plt.xticks(fontsize = 15) 
plt.yticks(fontsize = 15) 
plt.title("DR rates coefficents all configurations",fontsize=40)
plt.xlabel("electron energy (ev) ",fontsize=30)
plt.ylabel("cross section (barn)",fontsize=30)

plt.show()