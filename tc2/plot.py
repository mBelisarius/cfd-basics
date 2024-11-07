import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def solution_analytic(x, t, alpha, l):
    return np.sin(np.pi * x / l) * np.exp(-alpha * np.pi ** 2 * t / (l ** 2))


conditions = {'t': 20.0, 'alpha': 1.17e-4, 'l': 0.1}

df = pd.DataFrame({
    'P': [],
    'Xp': [],
    'TpAnalytic': [],
    'TpNumeric': [],
    'Error': [],
})

df.loc[0] = [0, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+0, np.nan]
df.loc[1] = [1, 5.00000000000000e-03, 1.55358419552192e-02, 1.519114550741685e-02, np.nan]
df.loc[2] = [2, 1.50000000000000e-02, 4.50867694105048e-02, 4.408642135705005e-02, np.nan]
df.loc[3] = [3, 2.50000000000000e-02, 7.02242897378960e-02, 6.866621111609576e-02, np.nan]
df.loc[4] = [4, 3.50000000000000e-02, 8.84877673043450e-02, 8.652447370547318e-02, np.nan]
df.loc[5] = [5, 4.50000000000000e-02, 9.80894456765170e-02, 9.591311795710193e-02, np.nan]
df.loc[6] = [6, 5.50000000000000e-02, 9.80894456765170e-02, 9.591311795710201e-02, np.nan]
df.loc[7] = [7, 6.50000000000000e-02, 8.84877673043450e-02, 8.652447370547318e-02, np.nan]
df.loc[8] = [8, 7.50000000000000e-02, 8.84877673043450e-02, 6.866621111609558e-02, np.nan]
df.loc[9] = [9, 8.50000000000000e-02, 4.50867694105048e-02, 4.408642135705026e-02, np.nan]
df.loc[10] = [10, 9.50000000000000e-02, 1.55358419552192e-02, 1.519114550741688e-02, np.nan]
df.loc[11] = [11, 1.00000000000000e-01, 0.00000000000000e+00, 0.00000000000000e+0, np.nan]

df['Error'] = df['TpAnalytic'] - df['TpNumeric']

xplot = np.linspace(0.0, conditions['l'], 100)

plt.scatter(df['Xp'], df['TpNumeric'], color='black', label='Simulation Data')
plt.plot(xplot, solution_analytic(xplot, *conditions.values()), color='red', label='Analytical Solution')

plt.grid()

# plt.title('Simulation Results')
plt.xlabel(r'Volume Centre $X_P$ $(m)$')
plt.ylabel(r'Temperature')
plt.legend()

plt.savefig('plot.png')
plt.show()
