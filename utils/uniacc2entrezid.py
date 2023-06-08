#!/usr/bin/env python

import click
import json
import re
import ssl
from urllib import parse, request

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "-u", "--uniacc", "uniaccs",
    help="UniProt Accession (e.g. P49711).",
    multiple=True,
    required=True,
)

def cli(**params):

    print("UniProt Accession\tEntrez Gene ID")
    for uniacc, entrezid in __get_uniaccs_entrezids(params["uniaccs"]):
        print("%s\t%s" % (uniacc, entrezid))

def __get_uniaccs_entrezids(uniaccs):


    return(uniaccs_entrezids)

def __get_uniaccs_entrezids(uniaccs):

    # Initialize
    uniaccs_entrezids = []
    gcontext = ssl.SSLContext()
    #url = "https://www.uniprot.org/uploadlists/"
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        #"from": "ACC+ID",
        "from": "UniProtKB_AC-ID",
        #"to": "P_ENTREZGENEID",
        "to": "GeneID",
        #"format": "tab",
        "query": " ".join(list(uniaccs))
    }

    # Get Entrez Gene IDs
    params = parse.urlencode(query=params).encode("utf-8")
    req = request.Request(url, params)
    with request.urlopen(req, context=gcontext) as handle:
        response = handle.read().decode("utf-8")
        response = json.loads(response)
        jobid = response["jobId"]
    url_jobid = f"https://rest.uniprot.org/idmapping/results/{jobid}"
    req_jobid = request.Request(url_jobid)
    with request.urlopen(req_jobid, context=gcontext) as handle:
        result = handle.read().decode("utf-8")
        result = json.loads(result)
    #for line in response.split("\n"):
    #    m = re.search("(\w+)\t(\d+)", line)
    #    if m:
    #        uniacc, entrezid = line.split("\t")
    #        uniaccs_entrezids.append([uniacc, int(entrezid)])
    for r in result["results"]:
        if "to" in r.keys():
            uniaccs_entrezids.append([r["from"], int(r["to"])])

    return(uniaccs_entrezids)

if __name__ == "__main__":
    cli()
