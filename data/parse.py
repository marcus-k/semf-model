import pandas as pd

# https://www-nds.iaea.org/amdc/
filename = "mass_1.mas20"
df = pd.read_fwf(
    filename,
    widths=(1,3,5,5,5,1,3,4,1,14,12,13,1,10,1,2,13,11,1,3,1,13,12),
    header=28, 
	index_col=False,
    usecols=(2, 3, 11, 13),
    names=("N", "Z", "E", "U(E)")
)

# Add useful columns from the dataset
df["known"] = 1 - pd.isnull(pd.to_numeric(df["E"], errors="coerce"))
df["E"] = pd.to_numeric(df["E"].str.replace("#",""))
df["U(E)"] = pd.to_numeric(df["U(E)"].str.replace("#",""))
df["A"] = df["N"] + df["Z"]
df["N-Z"] = df["N"] - df["Z"]

# If A if odd, subtract 1 from N-Z
# df.loc[df["A"] % 2 == 1, "N-Z"] -= 1

df["N/Z"] = df["N"] / df["Z"]

# Save binding energies to CSV
df.to_csv("binding_energy_per_A.csv", index=False)