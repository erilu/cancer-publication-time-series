
#Erick Lu
#Date Extractor

#This program will open a csv file full of abstracts and extract the date, and
#write the output into another csv file.

import csv
import re
import urllib
import os

#Open the file
masterfile = open ('full_abstracts_cancer.txt') 

data = masterfile.read()

splitdata = data.split("\n\n")

#print (len(splitdata))


datefile = open ('cancer_dates.txt', 'w')

for row in splitdata:
    dates = re.findall(  ".\s(\d{4})\s(\w{3})"  , row)
    try:
        #print (dates[0])
        datefile.write(",".join(dates[0]))
        datefile.write ("\n")
    except:
        print (row)

datefile.close()
    
masterfile.close()
