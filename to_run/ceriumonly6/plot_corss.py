import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


print('ede')
list_confi=['5p 5p','5p 6s','5p 4f','5p 5d','5p 6p','5p 7s','5p 5f','5p 6d','5p 7p', '5p 8s', '5p 5g']


convolution=[]
for i in range (1,len(list_confi)):

    y_1 = np.loadtxt('/home/flury/Fileserver/test_big_calculs/n=%s/arr_cross_strontum_%s.csv'%(i,i), delimiter=' ')
    

    x_1=y_1[:,0]
    y_1=y_1[:,1]
    convolution.append(y_1)

    ind=np.where(y_1< 1)[0]
    x_1=np.delete(x_1,ind)
    y_1=np.delete(y_1,ind)

    plt.plot(x_1, y_1,label="configuration %s"%(list_confi[i]))


'''    
print(convolution)
conv=convolution[0]
for i in range(1,len(convolution)):

    s=np.sum(convolution[i])
    conv=signal.fftconvolve(conv,convolution[i]/s,'same')#convolution is commutative


    conv=conv+convolution[i]

plt.plot(x_1, conv,label='sum_results')
'''





'''
x_2=y_2[:,0]
y_2=y_2[:,1]
ind=np.where(y_2< 0.01)[0]
x_2=np.delete(x_2,ind)
y_2=np.delete(y_2,ind)
#plt.plot(x_2, y_2,label="δn=6")


x_3=y_3[:,0]
y_3=y_3[:,1]
ind=np.where(y_3< 0.01)[0]
x_3=np.delete(x_3,ind)
y_3=np.delete(y_3,ind)
#plt.plot(x_3, y_3,label="δn=7")

x_4=y_4[:,0]
y_4=y_4[:,1]
ind=np.where(y_4< 0.01)[0]
x_4=np.delete(x_4,ind)
y_4=np.delete(y_4,ind)
#plt.plot(x_4, y_4,label="δn=8")

x_5=y_5[:,0]
y_5=y_5[:,1]
ind=np.where(y_5< 0.01)[0]
x_5=np.delete(x_5,ind)
y_5=np.delete(y_5,ind)
plt.plot(x_5, y_5,label="δn=9")


x_6=y_6[:,0]
y_6=y_6[:,1]
ind=np.where(y_6< 0.01)[0]
x_6=np.delete(x_6,ind)
y_6=np.delete(y_6,ind)
plt.plot(x_6, y_6,label="δn=10")


x_7=y_7[:,0]
y_7=y_7[:,1]
ind=np.where(y_7< 0.01)[0]
x_7=np.delete(x_7,ind)
y_7=np.delete(y_7,ind)
#plt.plot(x_7, y_7,label="δn=12")


'''
plt.legend(loc="upper right",fontsize=20)


plt.xticks(fontsize = 15) 
plt.yticks(fontsize = 15) 

plt.title("cross sections all configurations",fontsize=40)
plt.xlabel("electron energy (ev) ",fontsize=30)
plt.ylabel("cross section (barn)",fontsize=30)

plt.show()