#!/usr/bin/env bash

URL="https://stringdb-static.org/download/"

PREFIX="protein.links"
if ! [ -f "9606.${PREFIX}.v11.0.txt.gz" ]; then
    curl -L -O "${URL}/${PREFIX}.v11.0/9606.${PREFIX}.v11.0.txt.gz"
fi

PREFIX="protein.info"
if ! [ -f "9606.${PREFIX}.v11.0.txt.gz" ]; then
    curl -L -O "${URL}/${PREFIX}.v11.0/9606.${PREFIX}.v11.0.txt.gz"
fi