#!/usr/bin/env python3
# Generate a list of all characters that appear in
# the final (filtered) TF-IDF lists

import os
import gzip
import numpy as np

filesPath = "../dbTF/gene2wordclouds/filtered"

charactersArray = []

files = os.listdir(filesPath)

for fileName in files:
	
	filePath = os.path.join(filesPath, fileName)
	
	file = gzip.open(filePath, "rt")
	
	for line in file:
		
		line = line.rstrip("\n")
		line = line.split("\t")
		word = line[0]
		
		for character in word:
			
			charactersArray.append(character)
	
	file.close()



uniqueCharacters = np.unique(charactersArray)

outFile = open("characters_list.txt", "w")

for character in uniqueCharacters:
	outFile.write(character)

outFile.close()