import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import time

for i in range(len(sys.argv)):
    if (sys.argv[i] == "-i"): fileName = sys.argv[i+1]

startT = time.time()

data = np.loadtxt(fileName, dtype = int)
data = pd.DataFrame(data, columns = ['timestep', 'mol_id',
                                     'type_1', 'type_2', 'type_3'])
data['molecule'] = pd.Series([f"C{data.iloc[i,2]:>02}H{data.iloc[i,3]:>02}O{data.iloc[i,4]:>02}" for i in range(data.shape[0])])
data['time'] = data.timestep * 0.1 / 1000000

counts = data.groupby('molecule').time.value_counts().sort_index()

# plt.plot(counts['C01H04O00'], label = 'CH4')
# plt.plot(counts['C00H00O02'], label = 'O2')
plt.plot(counts['C01H03O00'], label = 'CH3', alpha = 0.3)
plt.plot(counts['C00H01O01'], label = 'OH' , alpha = 0.3)
plt.plot(counts['C00H01O02'], label = 'HOO', alpha = 0.3)
plt.plot(counts['C01H02O01'], label = 'CH2O',alpha = 0.3)
plt.plot(counts['C01H01O01'], label = 'CHO', alpha = 0.3)
plt.title("Time dependences of the numbers of main molecules/radicals")
plt.xlabel("Time (ns)")
plt.ylabel("Number of molecules")
plt.legend()
plt.grid()
plt.savefig(fileName[:-3]+"_molecules.png")
plt.show()

endT = time.time()
print(f"--- RUNTIME : {endT - startT} seconds ---")
