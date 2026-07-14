import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

data_a = pd.read_csv('Files/Data0MeV.csv', names=['centrality', 'dN_deta'])
data_c = pd.read_csv('Files/Data270MeV.csv', names=['centrality', 'dN_deta'])

x_new = np.linspace(0, 80, 40)
interp_a = interp1d(data_a['centrality'], data_a['dN_deta'], kind='linear', fill_value="extrapolate")
interp_c = interp1d(data_c['centrality'], data_c['dN_deta'], kind='linear', fill_value="extrapolate")

y_a_new = interp_a(x_new)
y_c_new = interp_c(x_new)


ratio = y_a_new / y_c_new


plt.figure(figsize=(7, 5))
plt.text(0.15, 0.9, 
         "ATLAS Internal\n$O+O$ 5.36 TeV", 
         transform=plt.gca().transAxes,  
         fontsize=14, 
         horizontalalignment='center',  
         verticalalignment='center')
plt.xlim(80, 0)
plt.plot(x_new, ratio, 'o-', color='red', label='O+O: (pT > 0) / (pT > 270)')
plt.xlabel('Centrality [%]')
plt.ylabel('Extrapolation')
plt.grid(True, linestyle=':', alpha=0.6)
plt.ylim(1.3, 1.8)
plt.legend()
plt.savefig("Plots/Extrapolation.png")