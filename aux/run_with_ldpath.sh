#!/bin/sh

here=`dirname $0`
LD_LIBRARY_PATH="$here"/libs:"$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH
exec "$0".bin "$@"
