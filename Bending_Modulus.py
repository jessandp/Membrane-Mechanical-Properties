import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

# Path to the PDB file
file_path = "/volumes/Macintosh HD/Users/jessicapallarez/Desktop/Sizes/Lipids/T1_POPC_200/gromacs/PO4_POPC_200.pdb"

# Read PDB file
d = pd.read_csv(file_path, sep=" ", comment="HEADER", skip_blank_lines=True, header=None)
d.columns = ["record", "atom", "atom_name", "residue_name", "chain_id", "residue_number", "x", "y", "z"]

##########################################################################################

# Identifying the lipid type
Lipid_Type = d["residue_name"].iloc[0]

if Lipid_Type == "DOPC":
    print("Lipid Type: DOPC")
else:
    print("Not DOPC!")

##########################################################################################

# Select coordinate values and create new DataFrame "d2" with selected columns
d2 = d[["x", "y", "z"]].copy()

# Count number of snapshots
n = len(d2) // 2

# Identifying upper leaflet and lower leaflet
leafs = np.arange(n) < n/2

# Determining cell length and height
cell_length = d.iloc[0]["x"]
cell_height = d.iloc[0]["z"]

# Cell length in nm
cell_length_nm = cell_length / 10

# Convert coordinate values to nm
d2["x"] = (d2["x"] % cell_length) / 10
d2["y"] = (d2["y"] % cell_length) / 10
d2["z"] = (d2["z"] % cell_height) / 10

# Separating upper and lower leaflet
P_split = np.split(d2, np.arange(n, len(d2), n))

# Summarizing the coordinate snapshots
summary = P_split[0].describe()

# Parameters for calculating q values
q_Nmax = 20
q_range = np.arange(-q_Nmax, q_Nmax) * 2 * math.pi / cell_length_nm
q_max = abs(q_range[0]) * math.sqrt(2)
q_0_05 = np.arange(0, q_max, 0.05)

# Creating DataFrame with q and v values
dd = pd.DataFrame({"qx": np.tile(q_range, len(q_range)), "qy": np.repeat(q_range, len(q_range))})
dd["q"] = np.sqrt(dd["qx"]**2 + dd["qy"]**2)
dd["v"] = np.searchsorted(q_0_05, dd["q"], side="right") - 1

##########################################################################################

# Fourier Transform
def fourier_transform(xyz):
    xyz["z"] = xyz["z"] - xyz["z"].mean()
    m = np.outer(dd["qx"], xyz["x"]) + np.outer(dd["qy"], xyz["y"])
    u = np.exp(-1j * m).dot(xyz["z"]) / (2 * n)
    u2 = np.real(np.conj(u) * u)
    return u2

# Perform Fourier Transform on upper leaflet coordinates
xyz = P_split[0].copy()
xyz["z"] = xyz["z"] - xyz["z"].mean()
m = np.outer(dd["qx"], xyz["x"]) + np.outer(dd["qy"], xyz["y"])
u = np.exp(-1j * m).dot(xyz["z"]) / (2 * n)
dd["u"] = u

# Apply Fourier Transform function to all leaflet snapshots
z = np.apply_along_axis(fourier_transform, axis=1, arr=P_split)

# Sort and find unique v values
unique_v = np.sort(np.unique(dd["v"]))

# Find indices of v = 1
indices_v1 = np.where(dd["v"] == 1)[0]

# Plotting ave_u2 (q) on a double log plot
fig, ax = plt.subplots()
ax.loglog(dd["q"][1:], dd["ave_u2"][1:], "o")
ax.set_xlim(left=1e-7, right=1)
plt.show()

##########################################################################################

# Finding the intercepts using the first 3 points
z3 = pd.DataFrame({"q": dd["q"], "u2": dd["ave_u2"]}).groupby("v").mean().reset_index()
z3["lg_q"] = np.log10(z3["q"])
z3["lg_u2"] = np.log10(z3["u2"])
intercepts = np.polyfit(z3["lg_q"][1:4], z3["lg_u2"][1:4], deg=1)

# First intercept (a)
a = intercepts[1]

# Second intercept (b)
b = intercepts[0]

##########################################################################################

# Plotting ave_u2 (q) on a double log plot and saving as PDF
fig, ax = plt.subplots()
ax.loglog(dd["q"][1:], dd["ave_u2"][1:], "o")
ax.plot(dd["q"], 10**(a + b * np.log10(dd["q"])), color="red")
ax.set_xlim(left=1e-7, right=1)
plt.savefig("/volumes/Macintosh HD/Users/jessicapallarez/Desktop/Lab Research/Jessica/P-pdbs/FSA_Results.pdf")
plt.close()

##########################################################################################

# Finding the bending modulus
k_B = 1.38e-23  # Boltzmann Constant
Temperature = 303.15  # Temperature membrane was simulated in
bending_modulus = k_B * Temperature / cell_length_nm**2 / 10**a

print("Lipid Type:", Lipid_Type)
print("Bending Modulus:", bending_modulus)

