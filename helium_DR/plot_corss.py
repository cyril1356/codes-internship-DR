import numpy as np
import matplotlib.pyplot as plt

y_1 = np.loadtxt('arr_cross_strontum_2.csv', delimiter=' ')
'''
y_2 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_2.csv', delimiter=' ')
y_3 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_3.csv', delimiter=' ')
y_4 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_4.csv', delimiter=' ')
y_5 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_5.csv', delimiter=' ')
y_6 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_6.csv', delimiter=' ')
y_7 = np.loadtxt('/home/flury/Fileserver/strontum5P/arr_cross_strontum_7.csv', delimiter=' ')
'''
x_1=y_1[:,0]
y_1=y_1[:,1]
x_1=x_1
ind=np.where(y_1< 0.01)[0]
x_1=np.delete(x_1,ind)
y_1=np.delete(y_1,ind)
plt.plot(x_1, y_1)


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


plt.title("cross sections helium",fontsize=40)
plt.xlabel("electron energy (ev) ",fontsize=30)
plt.ylabel("cross section (barn)",fontsize=30)
plt.savefig('all corss section ultime fig',dpi=30)
plt.show()