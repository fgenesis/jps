#!/bin/sh
echo "*** Waypoint path: ***"
./main
echo
echo "*** Detailed path: ***"
./main 1
echo
echo "*** Every-2nd-step path: ***"
./main 2
