import csv
import os

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python'
for filename in os.listdir(dir):
    data = []
    if filename.startswith("readshock") and filename.endswith(".csv"):
        # read in all the data
        with open(filename,'r') as csvFile:
            reader = csv.reader(csvFile)
            for row in reader:
                data.append(row)

        # make whatever changes you need
        data[0] = [800,250] + data[0]

        # now overwrite the old file
        with open(filename,'w',newline='') as csvFile:
            writer = csv.writer(csvFile)
            for row in data:
                writer.writerow(row)
            