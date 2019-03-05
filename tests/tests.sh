#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    [ "$MYTMP" = "" ] || rm -r "${MYTMP}"
    exit 1
}

[ $# -ne 6 ] && exit 1

SRCDIR=$1
BINDIR=$2
CASE=$3
INDEX_FLAGS=$4
FLAGS=$5
EXPECTED_FOLDER=$6
