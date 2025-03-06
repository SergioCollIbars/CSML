import re

def format_string(s):
    match = re.match(r'RECOEF\s*(\d{1,3})\s*(\d{1,3})(.*)', s)
    if match:
        formatted = f"RECOEF {match.group(1)} {match.group(2)}{match.group(3)}".rstrip() + '\n'
        return formatted
    return None  # Return None if the format doesn't match

fileName = 'GGM05S.txt'
fileName2 = "parsedFile.txt"

file = open(fileName, 'r')
lines = file.readlines()

file1 = open(fileName2, "a")  # append mode
file1.write("Today \n")

for line in lines:
    formatted_s = format_string(line)
    if formatted_s == None:
        s = line
    else: s = formatted_s
    file1.write(s)
file1.close()
file.close()
