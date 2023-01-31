# Import libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set working directory
os.chdir("output/data")

# Read the file stress_at_time0.dat
stress_at_time0 = pd.read_csv("stress_at_time0.dat", sep = "\t", header = None)

# rename Columns
stress_at_time0.columns = ["x", "y", "z", "s11", "s22", "s12", "s13", "s33", "s23"]

# filter the data for y = 0 and z = 0
stress_at_time0 = stress_at_time0[(stress_at_time0.y == 0) & (stress_at_time0.z == 0)]

# Plot the stress in x direction
plt.plot(stress_at_time0.x, stress_at_time0.s11)
plt.xlabel("x")
plt.ylabel("s11")
plt.show()