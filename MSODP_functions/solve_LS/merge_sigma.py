from uniplot import plot
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("file_1", help="The path to the data folder")
    parser.add_argument("file_2", help="plot title")
    
    args = parser.parse_args()

    file_1     = str(args.file_1)
    file_2     = str(args.file_2)
    
    matrix_1   = np.loadtxt(file_1)
    N_points   = len(matrix_1)
    data_1     = np.abs(matrix_1.reshape(N_points, 1))

    matrix_2   = np.loadtxt(file_2)
    N_points   = len(matrix_2)
    data_2     = np.abs(matrix_2.reshape(N_points, 1))

    y = data_1.T / data_2.T
    np.savetxt("data/sigma_ratio.txt", y, fmt='%.18e', delimiter='  ')

if __name__ == "__main__":
    main()