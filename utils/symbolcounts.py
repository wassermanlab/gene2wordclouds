import sys
import os
import numpy as np
import pandas as pd
import gzip
import re

def symbolcount():
    with gzip.open(sys.argv[1], 'rt') as file:
        global idfs
        idfs = pd.read_table(file)
        file.close

    cwd = os.getcwd()
    outdir = os.path.join(cwd, 'symbolstats')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    os.chdir(outdir)

    countstart(idfs)
    countend(idfs)
    summarize()

def countstart(input):
    #obtain words with non-alphabet & non-number start
    words_symbolstart = []
    for word in input["Word"]:
        if not re.findall("^[a-z0-9]", word):
            words_symbolstart.append(word)

    #obtain the starting symbol of the word
    startsymbol = []
    for word in words_symbolstart:
        startsymbol.append(word[0])

    #output file to show all words with symbol start and what their symbol is
    df_startsymbol = pd.DataFrame()
    df_startsymbol["Words"] = words_symbolstart
    df_startsymbol["Start"] = startsymbol
    df_startsymbol.to_csv('wordswithstartsymbol.txt', sep="\t", index=False)

    #get count of each symbol
    global df_start_sum
    df_start_sum = pd.DataFrame(df_startsymbol["Start"].value_counts())

def countend(input):
    words_symbolend = []
    for word in idfs["Word"]:
        if not re.findall("[a-z0-9]$", word):
            words_symbolend.append(word)

    endsymbol = []
    for word in words_symbolend:
        endsymbol.append(word[-1])

    #output file of words with ending symbol
    df_endsymbol = pd.DataFrame()
    df_endsymbol["Words"] = words_symbolend
    df_endsymbol["End"] = endsymbol
    df_endsymbol.to_csv('wordswithendsymbol.txt', sep="\t", index=False)

    global df_end_sum
    df_end_sum = pd.DataFrame(df_endsymbol["End"].value_counts())

def summarize():
    df_final_sum = df_start_sum.join(df_end_sum).fillna(0).astype(int)
    df_final_sum.to_csv('summarizesymbols.txt', sep="\t", index=True)

if __name__ == "__main__":
    symbolcount()
