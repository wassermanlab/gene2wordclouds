#!/usr/bin/env python

from Bio import Entrez, Medline
import click
import hashlib
import os

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "--dummy-dir",
    help="Dummy directory.",
    default="/tmp/",
    show_default=True,
)
@click.option(
    "-e", "--email",
    help="E-mail address.",
    default="Your.Name.Here@example.org",
    show_default=True,
)
@click.option(
    "-p", "--pmid", "pmids",
    help="PubMed ID (e.g. 8649389).",
    type=int,
    multiple=True,
    required=True,
)

def pmid2abstract(pmids, email, dummy_dir):

    print("PubMed ID\tAbstract")
    for pmid, abstract in __get_pmids_abstracts(pmids, email):
        print("%s\t%s" % (pmid, abstract))

def __get_pmids_abstracts(
    pmids, email="Your.Name.Here@example.org", dummy_dir="/tmp/"
):

    # Initialize
    ids = ",".join(list(map(str, sorted(pmids))))
    pmids_abstracts = {}
    for pmid in pmids:
        pmids_abstracts.setdefault(pmid, None)

    # Get abstracts
    h = hashlib.md5()
    h.update(ids.encode("utf-8"))
    prefix = h.hexdigest()
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)
    dummy_file = os.path.join(dummy_dir, "%s.%s.txt" % (
        os.path.basename(__file__), prefix)
    )
    Entrez.email = email
    handle = Entrez.efetch(
        db="pubmed", id=ids, rettype="medline", retmode="text"
    )
    data = handle.read()
    with open(dummy_file, "w") as handle:
        handle.write(data)
    with open(dummy_file) as handle:
        records = Medline.parse(handle)
        for record in records:
            if "PMID" in record:
                if "AB" in record:
                    pmids_abstracts[int(record["PMID"])] = record["AB"]
    os.remove(dummy_file)

    return([[k, v] for k, v in pmids_abstracts.items()])

if __name__ == "__main__":
    pmid2abstract()