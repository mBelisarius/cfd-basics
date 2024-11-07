import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def solution_analytic(t, alpha, l):
    return (2.0 / np.pi) * np.exp(-alpha * np.pi ** 2 * t / (l ** 2))


conditions = {'t': np.linspace(0.0, 20.0, 6), 'alpha': 1.17e-4, 'l': 0.1}

df = pd.DataFrame({
    'P': [],
    'Xp': [],
    'TpAnalytic': [],
    'TpNumeric': [],
    'Error': [],
})

df.loc[0] = [0, 0.00000000000000e+00, 6.36619772367581e-01, 6.31423598897955e-01, np.nan]
df.loc[1] = [1, 4.00000000000000e+00, 4.01125797554293e-01, 3.96070435305087e-01, np.nan]
df.loc[2] = [2, 8.00000000000000e+00, 2.52744122076472e-01, 2.48441442474679e-01, np.nan]
df.loc[3] = [3, 1.20000000000000e+01, 1.59250767798250e-01, 1.55838822686563e-01, np.nan]
df.loc[4] = [4, 1.60000000000000e+01, 1.00341827283559e-01, 9.77523653639599e-02, np.nan]
df.loc[5] = [5, 2.00000000000000e+01, 6.32240738415718e-02, 6.13167166532567e-02, np.nan]

df['Error'] = df['TpAnalytic'] - df['TpNumeric']

xplot = np.linspace(0.0, 20.0, 60)

plt.scatter(df['Xp'], df['TpNumeric'], color='black', label='Simulation Data')
plt.plot(xplot, solution_analytic(xplot, conditions['alpha'], conditions['l']), color='red',
         label='Analytical Solution')

plt.grid()

# plt.title('Simulation Results')
plt.xlabel(r'Time $t$ $(s)$')
plt.ylabel(r'Average Temperature')
plt.legend()

plt.savefig('plot2.png')
plt.show()
