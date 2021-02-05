#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import gzip
import matplotlib.pyplot as plt
import re
from wordcloud import WordCloud

# Compile
global regexp
regexp = re.compile("[^a-zA-Z]")

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@optgroup.group("Words", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--word", "words",
    help="Word (e.g. \"11-zinc-finger\").",
    multiple=True
)
@optgroup.option(
    "--words-file",
    help="File with words.",
    type=click.File("rt")
)
@optgroup.group("Weights", cls=RequiredMutuallyExclusiveOptionGroup)
@optgroup.option(
    "--weight", "weights",
    help="Weight (e.g. 1.0).",
    type=float,
    multiple=True
)
@optgroup.option(
    "--weights-file",
    help="File with weights.",
    type=click.File("rt")
)
@click.option(
    "-o", "--output-file",
    help="Output file.",
    type=click.Path(),
    required=True
)
@click.option(
    "--numwords",
    help="Number of words to display on word cloud.",
    type=int,
    default=200,
    show_default=True
)


def cli(**params):

    # Parse words
    if params["words_file"] is not None:
        params["words"] = []
        if params["words_file"].name.endswith(".gz"):
            handle = gzip.open(params["words_file"].name, "rt")
            params["words_file"].close()
        else:
            handle = params["words_file"]
        for line in handle:
            params["words"].append(line.strip("\n"))
        handle.close()

    # Parse weights
    if params["weights_file"] is not None:
        params["weights"] = []
        if params["weights_file"].name.endswith(".gz"):
            handle = gzip.open(params["weights_file"].name, "rt")
            params["weights_file"].close()
        else:
            handle = params["weights_file"]
        for line in handle:
            params["weights"].append(int(line.strip("\n")))
        handle.close()

    __get_word_cloud(params["words"], params["weights"], params["output_file"], params["numwords"])

def __get_word_cloud(words, weights, output_file, numwords):

    # Initialize
    frequencies = {}

    # Make word cloud
    for word, weight in zip(words, weights):
        frequencies.setdefault(word, 0.)
        frequencies[word] += weight
    wc = WordCloud(background_color="white", max_words=numwords)
    wc.generate_from_frequencies(frequencies)
    plt.imshow(wc, interpolation="bilinear")
    plt.axis("off")
    plt.savefig(output_file)

if __name__ == "__main__":
    cli()