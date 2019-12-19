#!/usr/bin/env bash

# For every argument, cat it into a temp file (which gets deleted later)
# then use head to pull the header off, print it, and grep to strip out
# every line with the same header in the concatenated file.
# Everything gets printed to standard output

tempfile=$(mktemp /tmp/strip_headers.XXXXXX)
exec 3>"$tempfile"
rm "$tempfile"

cat $@ > $tempfile && \
HEADER=$(head -n1 $tempfile) && \
echo ${HEADER} && \
grep -v "^${HEADER}" $tempfile