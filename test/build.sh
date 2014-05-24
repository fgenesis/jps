#!/bin/sh
g++ main.cpp -DNDEBUG -o main -O2 -pipe -Wall -pedantic
g++ main2.cpp ScenarioLoader.cpp -DNDEBUG -o main2 -O2 -pipe -Wall -pedantic
