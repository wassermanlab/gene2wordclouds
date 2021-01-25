#!/usr/bin/env python

import click
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

def uniacc2entrezid(uniaccs):

    print("UniProt Accession\tEntrez Gene ID")
    for uniacc, entrezid in __get_uniaccs_entrezids(uniaccs):
        print("%s\t%s" % (uniacc, entrezid))

def __get_uniaccs_entrezids(uniaccs):

    # Initialize
    uniaccs_entrezids = []
    gcontext = ssl.SSLContext()
    url = "https://www.uniprot.org/uploadlists/"
    params = {
        "from": "ACC+ID",
        "to": "P_ENTREZGENEID",
        "format": "tab",
        "query": " ".join(list(uniaccs))
    }

    # Get Entrez Gene IDs
    params = parse.urlencode(query=params).encode("utf-8")
    req = request.Request(url, params)
    with request.urlopen(req, context=gcontext) as handle:
        response = handle.read().decode("utf-8")
    for line in response.split("\n"):
        m = re.search("(\w+)\t(\d+)", line)
        if m:
            uniacc, entrezid = line.split("\t")
            uniaccs_entrezids.append([uniacc, int(entrezid)])

    return(uniaccs_entrezids)

if __name__ == "__main__":
    uniacc2entrezid()