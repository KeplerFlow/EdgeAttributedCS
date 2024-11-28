## Introduction
The problem of LCQ
- Code author: Liu Chuanyu from Hunan University

## File Structure Overview
- **Algorithm** Folder
  - **CoreGroup**: Core decomposition algorithm
  - **Online**: Online search algorithm
- **BuildIndex** Folder
  - **hierarchy**: Core list construction and search algorithm
  - **LCS**: Index algorithm, includes both basic and optimized index construction
- **core-maint** Folder
  - **glist** Folder: Order-based core maintenance algorithm
- **Dataset** Folder
  - Stores datasets required for experiments, including real datasets from Sem_network series and Movie series as well as other synthetic datasets

- **Graph.h/cpp**: Graph data structure definition
- **executor.cpp**: Label adjustment experiment and k-value adjustment experiment
- **radio.cpp**: Scalability experiment
- **run.cpp**: casestudy1 experiment
- **baseline.cpp**: casestudy2 experiment

## Datasets

| Dataset      | Node Count | Edge Count  | Kmax | Max Degree | Nature          |
|--------------|------------|-------------|------|------------|-----------------|
| orkut        | 3,072,441  | 117,185,083 | 253  | 33,313     | Social Network  |
| twitter      | 456,626    | 14,855,842  | 125  | 51,386     | Social Network  |
| wikitalk     | 1,140,149  | 7,833,140   | 124  | 141,951    | Wikipedia talk  |
| com-youtube  | 1,134,890  | 2,987,624   | 51   | 28,754     | Video           |
| dblp.ungraph | 317,080    | 1,049,866   | 113  | 343        | Academic Papers |
| Sem_10M      | 2,135,390  | 10,080,227  | 54   | 81         | Academic Papers |
| Sem_2M       | 896,415    | 2,336,696   | 68   | 93         | Academic Papers |
| Sem_10W      | 48,131     | 206,950     | 76   | 98         | Academic Papers |
| Movie_2M     | 149,596    | 2,689,454   | 30   | 87         | Movies          |
| Movie_10W    | 597        | 113,658     | 201  | 515        | Movies          |

## Compilation Method
- `$ make`
  - Simply run the make command, to conduct different experiments just modify the first cpp file in the makefile.

## Execution Parameters
- The first parameter is the executable file name, determined by makefile.
- The second parameter is the dataset path.
- The third parameter is a boolean determining whether to use the optimized index construction method.
- The fourth parameter is a boolean determining whether to load an existing serialized index file.

## Execution Command
```
./a.out ./Dataset/<the dataset you need> <the bool of whether use our index> <the bool of whether use the serialized index file>

For example:
./a.out ./Dataset/orkut.txt true true
```

## Acknowledgments
Many thanks to seniors like Dong Wen, Ying Zhang, Yikai Zhang, Francesco Gullo for their support to the author.

## Contact
If you have any questions about the code, please email keplerliu10@gmail.com or 1783949591@hnu.edu.cn for consultation.
