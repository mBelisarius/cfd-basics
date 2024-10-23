import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def solution_analytic(x, ta, tb, k, qdot, l):
    return ta + x * (tb - ta) / l + x * (l - x) * qdot / (2.0 * k)


conditions = {'ta': 0.0, 'tb': 100.0, 'k': 400.0, 'qdot': 5.0e5, 'l': 1.0}

df = pd.DataFrame({
    'P': [],
    'Xp': [],
    'TpAnalytic': [],
    'TpNumeric': [],
    'Error': [],
})

df.loc[0] = [0, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00]
df.loc[1] = [1, 5.00000000000000e-02, 3.46875000000000e+01, 3.62500000000000e+01, -1.56249999999999e+00]
df.loc[2] = [2, 1.50000000000000e-01, 9.46875000000000e+01, 9.62500000000000e+01, -1.56250000000001e+00]
df.loc[3] = [3, 2.50000000000000e-01, 1.42187500000000e+02, 1.43750000000000e+02, -1.56250000000000e+00]
df.loc[4] = [4, 3.50000000000000e-01, 1.77187500000000e+02, 1.78750000000000e+02, -1.56250000000003e+00]
df.loc[5] = [5, 4.50000000000000e-01, 1.99687500000000e+02, 2.01250000000000e+02, -1.56250000000000e+00]
df.loc[6] = [6, 5.50000000000000e-01, 2.09687500000000e+02, 2.11250000000000e+02, -1.56250000000003e+00]
df.loc[7] = [7, 6.50000000000000e-01, 2.07187500000000e+02, 2.08750000000000e+02, -1.56250000000000e+00]
df.loc[8] = [8, 7.50000000000000e-01, 1.92187500000000e+02, 1.93750000000000e+02, -1.56250000000000e+00]
df.loc[9] = [9, 8.50000000000000e-01, 1.64687500000000e+02, 1.66250000000000e+02, -1.56249999999994e+00]
df.loc[10] = [10, 9.50000000000000e-01, 1.24687500000000e+02, 1.26250000000000e+02, -1.56249999999993e+00]
df.loc[11] = [11, 1.00000000000000e+00, 1.00000000000000e+02, 1.00000000000000e+02, 0.00000000000000e+00]

xplot = np.linspace(0.0, conditions['l'], 100)

plt.scatter(df['Xp'], df['TpNumeric'], color='black', label='Simulation Data')
plt.plot(xplot, solution_analytic(xplot, *conditions.values()), color='red', label='Analytical Solution')

plt.grid()

# plt.title('Simulation Results')
plt.xlabel(r'Volume Centre $X_P$ $(m)$')
plt.ylabel(r'Temperature $(°C)$')
plt.legend()

plt.savefig('plot.png')
plt.show()
