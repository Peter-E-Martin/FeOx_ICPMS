# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:39:01 2022

@author: Peter
"""

import pandas as pd
from uncertainties import ufloat
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Uatnorm = (0.000000001/238)*6.022e23
Uratnorm = 139.52
Thatnorm = (0.000000001/232)*6.022e23
Thratnorm = 1e99
Smatnorm = (0.00000001/150.36)*(6.022E+23)

def exponential(x, a, b):
    return a*np.exp(b*x)

def make_ref(measured, iso, spikeiso):
    return ufloat(np.average((measured[iso]/measured[spikeiso])),
                  np.std((measured[iso]/measured[spikeiso]), ddof=1))

def make_ratio(samples, iso, iso_s, spikeiso, spikeiso_s):
    return [(ufloat(samples[iso].iloc[n], samples[iso_s].iloc[n])/
             ufloat(samples[spikeiso].iloc[n], samples[spikeiso_s].iloc[n]))
            for n in range(len(samples['Sample Name']))]

def ICPMS_mols(atnorm, ratnorm, norm, spkblk, Nb_blk, R_meas):
    atoms_shot = (atnorm*(ratnorm-norm))/(ratnorm*(norm-spkblk))
    Atoms = (spkblk-R_meas)/(R_meas/ratnorm-1)*atoms_shot
    Nbat_blk = (spkblk-Nb_blk)/(Nb_blk/ratnorm-1)*atoms_shot
    mols = (Atoms-Nbat_blk)/6.022e23
    return mols


def Sm_ICPMS_mols(Smatnorm, NdSm, SmNdRat, Nb_Smblk):
    Nd145_shot = (1/NdSm)*Smatnorm*0.1502
    atoms147 = SmNdRat*Nd145_shot
    Nbat_Smblk = Nd145_shot*Nb_Smblk
    mols_147Sm = (atoms147-Nbat_Smblk)/6.022e23
    return mols_147Sm

# fn = r'C:\Users\pemar\Documents\FlowersResearch\Sparrow\Sparrow-CU-TRaIL\TRaIL-Data\icpmsData\Example_directory\FULL.xlsx'
# fn = r'C:\Users\pemar\Documents\FlowersResearch\Hematite_DoubleDating\Analyses\2022_03_17_ICPMS\full.xlsx'
fn = r'C:\Users\pemar\Documents\FlowersResearch\Hematite_DoubleDating\Analyses\2022_03_03_ICPMS\full.xlsx'

# import data
data = pd.read_excel(fn,
                     header = None)

data.iloc[0:2] = data.iloc[0:2].fillna(method='ffill', axis = 1)
data.iloc[0:2] = data.iloc[0:2].fillna('')
data.columns = data.iloc[0:2].apply(lambda x: '.'.join([y for y in x if y]), axis=0)
data = data.iloc[2:]

col_names = list(data.columns)
subs = {'Sample.Sample Name': 'Sample Name',
        '56  Fe  ': '56Fe', '57  Fe  ': '57Fe',
        '145  Nd  ': 'Nd145', '147  Sm  ': 'Sm147', '230  Th  ': 'Th230',
        '232  Th  ': 'Th232', '235  U  ': 'U235', '238  U  ': 'U238',
        '[ No Gas ]': '_ng', '[ He ]': '_he',
        ' .CPS': '', ' RSD': '_s', ' .Det.': '_det'}
for c in range(len(col_names)):
    for s in subs:
        if s in col_names[c]:
            col_names[c] = col_names[c].replace(s, subs[s])

to_replace = {k: v for k,v in zip(list(data.columns), col_names)}

data = data.rename(columns=to_replace)

# # separate blanks, normals, and samples
Nb_blanks = data.loc[data['Sample Name'].str.contains('~B')]
acid_blanks = data.loc[data['Sample Name'].str.contains('~A')]
norms = data.loc[data['Sample Name'].str.contains('~No')]
spikes = data.loc[data['Sample Name'].str.contains('~B')]

# # exclude all special analyses (incl. ~W which is the wash)
samples = data.loc[~data['Sample Name'].str.contains('~B') &
                    ~data['Sample Name'].str.contains('~A') &
                    ~data['Sample Name'].str.contains('~No') &
                    ~data['Sample Name'].str.contains('~S')&
                    ~data['Sample Name'].str.contains('~W')&
                    ~data['Sample Name'].str.contains('ppm')&
                    ~data['Sample Name'].str.contains('Water')&
                    ~data['Sample Name'].str.contains('Wash')]

# Regress Fe stds
iron_stds = data.loc[data['Sample Name'].str.contains('ppm')]

fe_iso = '56Fe_ng'

iron_stds.loc[:,'conc'] = iron_stds['Sample Name'].str.split('ppm').str[0].astype(int)
popt, pcov = curve_fit(exponential, iron_stds['conc'], iron_stds[fe_iso], p0=[1e7, 0.01])

x = np.linspace(0,500,1000)
y = exponential(x, *popt)

plt.scatter(iron_stds['conc'], iron_stds[fe_iso])
plt.plot(x, y)
plt.show()

# Get reduced data for U
# Get isotope ratio for blanks and normals
U_SpkBlnk = 0.00000001#make_ref(spikes, 'U238_ng', 'U235_ng')
U_NbBlnk = make_ref(Nb_blanks, 'U238_ng', 'U235_ng')
U_norm = make_ref(norms, 'U238_ng', 'U235_ng')

# Generate sample isotope ratios with uncertainty
samples.loc[:,'238R'] = make_ratio(samples, 'U238_ng', 'U238_ng_s',
                             'U235_ng', 'U235_ng_s')

# Calculate number of mols
mol238 = ICPMS_mols(Uatnorm, Uratnorm, U_norm, U_SpkBlnk, U_NbBlnk, samples['238R'])
samples.loc[:,'mol 238'] = [m.n for m in mol238]
samples.loc[:,'mol 238 s'] = [m.s for m in mol238]


# Repeat for Th
Th_SpkBlnk = make_ref(spikes, 'Th232_ng', 'Th230_ng')
Th_NbBlnk = make_ref(Nb_blanks, 'Th232_ng', 'Th230_ng')
Th_norm = make_ref(norms, 'Th232_ng', 'Th230_ng')

samples.loc[:,'232R'] = make_ratio(samples, 'Th232_ng', 'Th232_ng_s',
                             'Th230_ng', 'Th230_ng_s')

mol232 = ICPMS_mols(Thatnorm, Thratnorm, Th_norm, Th_SpkBlnk, Th_NbBlnk, samples['232R'])
samples.loc[:,'mol 232'] = [m.n for m in mol232]
samples.loc[:,'mol 232 s'] = [m.s for m in mol232]

# Repeat for Sm
# Get isotope ratio for blanks and normals
#TODO verify this with an apatite run
Sm_NbBlnk = make_ref(Nb_blanks, 'Sm147_ng', 'Nd145_ng')
Sm_norm = make_ref(norms, 'Sm147_ng', 'Nd145_ng')

# Generate sample isotope ratios with uncertainty
samples.loc[:,'147R'] = make_ratio(samples, 'Sm147_ng', 'Sm147_ng_s',
                              'Nd145_ng', 'Nd145_ng_s')

# Calculate number of mols
mol147 = Sm_ICPMS_mols(Smatnorm, Sm_norm, samples['147R'], Sm_NbBlnk)
samples.loc[:,'mol 147'] = [m.n for m in mol232]
samples.loc[:,'mol 147 s'] = [m.s for m in mol232]


# Calculate Fe std based concentrations
samples.loc[:,'ppm Fe'] = np.log((samples[fe_iso]/popt[0]).astype(float))/popt[1]
samples.loc[:,'g hematite'] = samples['ppm Fe']*(3/(0.7*1000000))

relevant_cols = ['Sample Name', 'mol 238', 'mol 238 s', 'mol 232', 'mol 232 s', 'mol 147', 'mol 147 s']
samples_trimmed = samples[relevant_cols]

# samples_trimmed.to_excel(r'C:\Users\pemar\Documents\FlowersResearch\Hematite_DoubleDating\Analyses\2022_03_17_ICPMS\init_results.xlsx')