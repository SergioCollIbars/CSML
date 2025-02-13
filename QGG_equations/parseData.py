
fileName = 'GRGM_400A.txt'
fileName2 = "parsedFile.txt"

file = open(fileName, 'r')
lines = file.readlines()

file1 = open(fileName2, "a")  # append mode
file1.write("Today \n")

for line in lines:
    nline = line[0:-3]
    nline = nline + '\n'
    file1.write(nline)
file1.close()
file.close()
