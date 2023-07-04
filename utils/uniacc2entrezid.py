#!/usr/bin/env python

import click
import time
from unipressed import IdMappingClient

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

    # Initialize
    uniaccs_entrezids = []

    # Get Entrez Gene IDs
    request = IdMappingClient.submit(
        source="UniProtKB_AC-ID", dest="GeneID", ids=set(uniaccs)
    )
    time.sleep(10) # wait for a few seconds
    results = list(request.each_result())
    for r in results:
        uniaccs_entrezids.append([r["from"], int(r["to"])])

    return(uniaccs_entrezids)

# def __get_uniaccs_entrezids(uniaccs):

#     # Initialize
#     uniaccs_entrezids = []
#     gcontext = ssl.SSLContext()
#     

#     # Do it in chunks of 25 accessions
#     for i in range(0, len(uniaccs), 25): 

#         # Initialize
#         query = {
#             "from": "UniProtKB_AC-ID",
#             "to": "GeneID",
#             "ids": " ".join(list(uniaccs[i:i+25]))
#         }

#         # Get Entrez Gene IDs
#         params = parse.urlencode(query=query).encode("utf-8")
#         req = request.Request(url, params)
#         with request.urlopen(req, context=gcontext) as handle:
#             response = handle.read().decode("utf-8")
#             response = json.loads(response)
#             jobid = response["jobId"]
#         url_jobid = f"https://rest.uniprot.org/idmapping/results/{jobid}"
#         req_jobid = request.Request(url_jobid)
#         with request.urlopen(req_jobid, context=gcontext) as handle:
#             result = handle.read().decode("utf-8")
#             result = json.loads(result)
#         for r in result["results"]:
#             if "to" in r.keys():
#                 uniaccs_entrezids.append([r["from"], int(r["to"])])

#         # Sleep between queries
#         time.sleep(1.)

#     return(uniaccs_entrezids)

if __name__ == "__main__":
    cli()
