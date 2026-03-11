import numpy as np

def load_matrix(filepath):

    vals = []

    with open(filepath, "r") as f:
        for line in f:

            if line.startswith("RECOEF"):

                # Remove keyword
                line = line.replace("RECOEF", "", 1)

                parts = line.split()

                # parts now: [degree, order, col4, col5, ...]
                if len(parts) >= 4:
                    c4 = float(parts[2].replace("D", "E"))
                    c5 = float(parts[3].replace("D", "E"))
                    vals.append([c4, c5])

    return np.array(vals)
