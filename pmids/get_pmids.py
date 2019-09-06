#This file, for each TAXON in JASPAR (i.e. vertebrates, insects, namatodes, fungi, plants), it fetches all TFs, then transforms this TFs to ENTREZ IDs, then identifies the PMIDs associated to individual ENTREZ IDs, then call the R script to get the word counts.

#For this to work, you need to:
#1) iterate the taxons.
#2) create a folder of the taxons.
#3) iterate through the individual TFs of each taxon.
#4) execute the R script for each PMID associated to each TF.

#Files that you need to download and store as JSONs for repeated use:
#Since JASPAR stores UniProt IDs, a UniProt ID to Entrez ID mapping.
#Since PMIDs are not going to change, download from NCBI the Entrez ID to PMID mapping.

import csv
import coreapi
import pickle
import os
import json

codec = coreapi.codecs.CoreJSONCodec()
client = coreapi.Client()

taxons = ["vertebrates", "plants", "insects", "nematodes", "fungi"]

for taxon in taxons:
    pubmed_ids = []
    try:
        file_name = taxon + "_uniprot.csv"
        open(file_name)

    except IOError as e:
        print(os.strerror(e.errno))
        file_name = taxon+"_uniprot.csv"
        csv_file = open(file_name, 'w')
        response = client.get("http://jaspar.genereg.net/api/v1/taxon/%s/" % taxon)
        json_obj = json.loads(codec.encode(response))
        while json_obj["next"] is not None:
            for profile in json_obj["results"]:
                if profile["collection"] == "CORE":
                    response2 = client.get("http://jaspar.genereg.net/api/v1/matrix/%s/" % profile["matrix_id"])
                    json_obj2 = json.loads(codec.encode(response2))
                    for uniprot in json_obj2["uniprot_ids"]:
                        for species in json_obj2["species"]:
                            csv_file.write(uniprot)
                            csv_file.write('\n')
            response = client.get(json_obj["next"])
            json_obj = json.loads(codec.encode(response))

        # Do last page
        for profile in json_obj["results"]:
            if profile["collection"] == "CORE":
                response2 = client.get("http://jaspar.genereg.net/api/v1/matrix/%s/" % profile["matrix_id"])
                json_obj2 = json.loads(codec.encode(response2))
                for uniprot in json_obj2["uniprot_ids"]:
                    for species in json_obj2["species"]:
                        csv_file.write(uniprot)
                        csv_file.write('\n')

    try:
        entrez_pubmed = pickle.load(open('gene2pubmed.p', 'rb'))

    except IOError as e:
        print(os.strerror(e.errno))
        linkage_dict = {}
        with open("gene2pubmed") as links:
            link_reader = csv.reader(links, delimiter='\t')
            for rownum, link in enumerate(link_reader):
                if rownum > 0:
                    temp = []
                    gene_id = link[1]
                    pm_id = link[2]
                    if gene_id in linkage_dict:
                        linkage_dict[gene_id].append(pm_id)
                    else:
                        temp.append(pm_id)
                        linkage_dict[gene_id] = temp
            pickle.dump(linkage_dict, open("gene2pubmed.p", "wb"), protocol=2)


    try:
        uniprot_entrez = pickle.load(open('uniprot_entrez.p', 'rb'))

    except IOError as e:
        print(os.strerror(e.errno))
        uniprot_entrez = {}
        with open("uniprot_to_entrez.csv") as links:
            link_reader = csv.reader(links, delimiter=',')
            for rownum, link in enumerate(link_reader):
                if rownum > 0:
                    temp = []
                    entrez_id = link[1]
                    uniprot_id = link[0]
                    if uniprot_id in uniprot_entrez:
                        uniprot_entrez[uniprot_id].append(entrez_id)
                    else:
                        temp.append(entrez_id)
                        uniprot_entrez[uniprot_id] = temp
            pickle.dump(uniprot_entrez, open("uniprot_to_entrez.p", "wb"), protocol=2)

    os.mkdir(taxon)
    with open(file_name) as uniprot_ids:
        enh_reader = csv.reader(uniprot_ids, delimiter='\t')
        for rownum, line in enumerate(enh_reader):
            uniprot_id=line[0]
            if uniprot_id in uniprot_entrez:
                entrez_id = uniprot_entrez[uniprot_id]
                entrez_id=entrez_id[0]
                if entrez_id in entrez_pubmed:
                    pmids= entrez_pubmed[entrez_id]
                    pubmed_ids = pubmed_ids + pmids
    pmid_set=set(pubmed_ids)
