#!/bin/bash

vulcanFlaresDir='vulcan-flares'
currentDir="$(basename "$PWD")"
if [ "$currentDir" != "$vulcanFlaresDir" ]; then
    echo "[ERROR] This script is meant to be run from $vulcanFlaresDir" >&2
    exit 1
fi

export VULCAN_RUNS_ARE_CHAINED=1

# if n=1 and n<3, then it means 2 times (3 - 1 = 2)
for ((n=1; n<3; n++)); do
    export VULCAN_RUN_ORDINAL_NUMBER=${n}
    echo -e "\n[INFO] VULCAN run #${VULCAN_RUN_ORDINAL_NUMBER}\n"
    python ./vulcan.py
    if [ $? -eq 0 ]; then
        echo -e "\n[INFO] Done executing VULCAN run #${VULCAN_RUN_ORDINAL_NUMBER}\n"
    else
        echo "[ERROR] The VULCAN run #${VULCAN_RUN_ORDINAL_NUMBER} failed" >&2
        exit 2
    fi
done
