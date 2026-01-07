import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# all csv paths after script name
files = sys.argv[1:]

if len(files) == 0:
    raise ValueError("Please provide at least one CSV file.")

plt.figure()

for fin in files:
    df = pd.read_csv(fin)
    raw = df[df["Type"] == "Raw"]

    label = os.path.basename(fin)

    plt.plot(raw["Volume"], raw["Energy"], "o-", label=label)

plt.xlabel("Volume")
plt.ylabel("Energy")
plt.legend()
plt.tight_layout()
plt.show()
