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
plt.xlim(0, 80)
plt.plot(data_a['centrality'],data_a['dN_deta'],  'o-', color='red', label='O+O:  (pT > 0)')
plt.plot(data_c['centrality'],data_c['dN_deta'],  'o-', color='green', label='O+O: (pT > 270) ')

plt.xlabel('Centrality [%]')
plt.ylabel('dn/deta')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig("Plots/ExtrapolationControl.png")


plt.figure(figsize=(7, 5))
plt.text(0.15, 0.9, 
         "ATLAS Internal\n$O+O$ 5.36 TeV", 
         transform=plt.gca().transAxes,  
         fontsize=14, 
         horizontalalignment='center',  
         verticalalignment='center')
plt.xlim(0, 80)
plt.plot(x_new, ratio, 'o-', color='red', label='O+O: (pT > 0) / (pT > 270)')
plt.xlabel('Centrality [%]')
plt.ylabel('Extrapolation')
plt.grid(True, linestyle=':', alpha=0.6)
plt.ylim(1.3, 1.8)
plt.legend()
plt.savefig("Plots/Extrapolation.png")

x_targets = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 
    7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5])
interp_ratio = interp1d(x_new, ratio, kind='linear', fill_value="extrapolate")
ratio_targets = interp_ratio(x_targets)

# 4. Wypisujemy wyniki w czytelny sposób (zaokrąglone do 3 miejsc po przecinku)
print("\n=== WYNIKI DLA TWOICH PUNKTÓW CENTRALNOŚCI ===")
print("  x (Cent)  |  ratio (Extrapolation)")
print("------------+------------------------")
for x, r in zip(x_targets, ratio_targets):
    print(f"   {x:<8} |  {r:.3f}")
print("==============================================")