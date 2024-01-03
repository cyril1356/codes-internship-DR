import numpy as np
import matplotlib.pyplot as plt


y = np.loadtxt('DR_rates_strontum_2.csv', delimiter=',')
x=np.linspace(20 ,10**5,500000)

fig = plt.figure()
ax = fig.add_subplot()
plt.plot(x, y)
 
ax.set_aspect('auto', adjustable='box')

plt.yscale('log')
plt.xscale('log')
plt.xlim(0, 10**5.5)

plt.xticks(fontsize = 15) 
plt.yticks(fontsize = 15) 

plt.title("DR rates coefficients of helium-like carbon ",fontsize=35)
plt.xlabel("temperature K",fontsize=22)
plt.ylabel("DR rate coefficients cm³/s",fontsize=22)

plt.savefig('DR_rates')
plt.show()