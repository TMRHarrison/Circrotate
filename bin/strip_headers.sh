#!/usr/bin/env bash

# For every argument, cat it into a temp file (which gets deleted later)
# then use head to pull the header off, print it, and grep to strip out
# every line with the same header in the concatenated file.
# Everything gets printed to standard output

# set up temp file
tempfile=$(mktemp /tmp/strip_headers.XXXXXX)
function finish {
    [ -e $tempfile ] && rm $tempfile
    exit
}
# on some error codes, delete the temp file
trap finish EXIT

cat $@ > $tempfile && \
HEADER=$(head -n1 $tempfile) && \
echo ${HEADER} && \
grep -v "^${HEADER}" $tempfile