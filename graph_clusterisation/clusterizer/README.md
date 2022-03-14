# hgmc
## About 
The application finds the partition of graph's vertices 
that maximizes specific metric function. You can take a look  
at it here: https://codeforces.com/contest/1378/problem/A1

## Requirements
1) c++ 14 compiler
2) Boost
3) OpenMP
4) CMake

## How to build
1) Clone the repo
2) Move to cloned directory
3) Create a folder named **build**
4) Move to **build** folder and run from the terminal:  
   cmake ..  
   make

## How to run
From the **build** folder:
1) Run from the terminal (for example):  
   ./app/run_app --input ../data/graph.in --output output.out --algo louvain
2) To see all possible arguments run:  
   ./app/run_app --help
