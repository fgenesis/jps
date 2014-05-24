#!/bin/sh

for f in maps/*.scen; do
    echo "$f"
    time (./main2 "$f" > /dev/null)
done

