import os
import re
import numpy as np
from scipy.io import savemat

folder = r"C:\Users\edric\Downloads\maps\maps"   # CHANGE THIS

# grab files named like map.4, map.20, map.296, etc.
files = []
for name in os.listdir(folder):
    m = re.match(r"^map\.(\d+)$", name)
    if m:
        timestep = int(m.group(1))
        files.append((timestep, name))

files.sort(key=lambda x: x[0])

timesteps = np.array([t for t, _ in files])

# read first file
first_path = os.path.join(folder, files[0][1])
A = np.loadtxt(first_path)

lat = A[:, 0]
lon = A[:, 1]

n_points = A.shape[0]
n_times = len(files)

X = np.zeros((n_points, n_times))
X[:, 0] = A[:, 2]

# read remaining files
for k, (timestep, name) in enumerate(files[1:], start=1):
    path = os.path.join(folder, name)
    A = np.loadtxt(path)

    # optional safety check
    if A.shape[0] != n_points:
        raise ValueError(f"{name} has different number of rows")

    X[:, k] = A[:, 2]

print("X shape:", X.shape)
print("lat shape:", lat.shape)
print("lon shape:", lon.shape)
print("timesteps:", timesteps)

output_file = os.path.join(folder, "DMD_Data.mat")

savemat(output_file, {
    "X": X,
    "lat": lat,
    "lon": lon,
    "timesteps": timesteps
})

print(f"Saved to: {output_file}")