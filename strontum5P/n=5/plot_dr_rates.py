import numpy as np
import matplotlib.pyplot as plt


y = np.loadtxt('/media/flury/ESD-USB/strontum5P/n=5/DR_rates_strontum_2.csv', delimiter=',')
x=np.linspace(1 ,10**3,50000)

plt.plot(x, y)

plt.yscale('log')
plt.xscale('log')
plt.title("DR rates")
plt.xlabel("temp(ev)")
plt.ylabel("DR rate coef")
plt.savefig('DR_rates')
plt.show()