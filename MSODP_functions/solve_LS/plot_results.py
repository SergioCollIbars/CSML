from uniplot import plot
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("data_file", help="The path to the data folder")
    parser.add_argument("title", help="plot title")
    
    args = parser.parse_args()

    file     = str(args.data_file)
    tt       = str(args.title)
    matrix   = np.loadtxt(file)
    N_points = len(matrix)
    data = np.abs(matrix.reshape(N_points, 1))

    y = data.T
    plot(y, y_as_log=True, title=tt)

if __name__ == "__main__":
    main()