#!/bin/sh
g++ main.cpp -o main -O2 -pipe -Wall
g++ main2.cpp ScenarioLoader.cpp -o main2 -O2 -pipe -Wall
