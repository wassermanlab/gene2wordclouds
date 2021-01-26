#!/usr/bin/env python

import click
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import collections
import gzip
import matplotlib.pyplot as plt
import nltk
import numpy as np
import os
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

def words2cloud(words, words_file, weights, weights_file, output_file):

    # Parse words
    if words_file is not None:
        words = []
        if words_file.name.endswith(".gz"):
            handle = gzip.open(words_file.name, "rt")
            words_file.close()
        else:
            handle = words_file
        for line in handle:
            words.append(line.strip("\n"))
        handle.close()

    # Parse weights
    if weights_file is not None:
        weights = []
        if weights_file.name.endswith(".gz"):
            handle = gzip.open(weights_file.name, "rt")
            weights_file.close()
        else:
            handle = weights_file
        for line in weights_file:
            weights.append(float(line.strip("\n")))
        handle.close()

    __make_word_cloud(words, weights, output_file)

def __make_word_cloud(words, weights, output_file):

    # Initialize
    frequencies = {}

    # Make word cloud
    for word, weight in zip(words, weights):
        frequencies.setdefault(word, 0.)
        frequencies[word] += weight
    wc = WordCloud(background_color="white")
    wc.generate_from_frequencies(frequencies)
    plt.imshow(wc, interpolation="bilinear")
    plt.axis("off")
    plt.savefig(output_file)

if __name__ == "__main__":
    words2cloud()