#!/bin/sh
c++ main.cpp -DNDEBUG -o main -O2 -pipe -Wall -pedantic
c++ main2.cpp ScenarioLoader.cpp -DNDEBUG -o main2 -O2 -pipe -Wall -pedantic
