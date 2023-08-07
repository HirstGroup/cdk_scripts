#!/bin/bash
# test command line arguments

set -e

scripts=$scripts

# default values
complex=complex
part=""
time=time
overwrite=NO

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -c|--complex)
    complex="$2"
    shift
    ;;
    -p|--part)
    part="$2"
    shift
    ;;
    -t|--time)
    time="$2"
    shift
    ;;
    -O)
    overwrite=YES
    ;;
    *)
    echo "Unknown argument: $1"
    exit 1
    ;;
esac
shift
done

echo "complex =" $complex "time =" $time "part =" \"$part\" "overwrite =" $overwrite