# Modular DAG sampling

The code used in the experiments of the UAI 2019 paper "Exact Sampling of Directed Acyclic Graphs from Modular Distributions" (Topi Talvitie, Aleksis Vuoksenmaa, Mikko Koivisto). The code was written by Aleksis Vuoksenmaa, and adapted for distribution by Topi Talvitie.

This is a tool for sampling random directed acyclic graphs (or DAGs) from a modular probability distribution, i.e. a distribution where the probability weight of a DAG is given by a product of per-node weights, each of which only depends on the set of parents for the node. As a special case, one can use this tool to sample DAGs uniformly at random from the space of all DAGs with n nodes. 

## Compile

In the root folder, run
```
make
```

## Usage

To see the description of the command line arguments, run the program without arguments

```
./sampler
```

### Uniform case

To sample 10 DAGs with 5 nodes from the uniform distribution, run

```
./sampler symmetric uniform 5 10
```

The program writes the sampled DAGs into the standard output stream and other information to standard error stream.

### Symmetric case

To sample DAGs with symmetric parent set weights, you first need a file that contains the weights for each parent set size. The file should first contain the desired number of nodes *n*, and after that for each parent set size from 0 up to *n*-1 specify the natural logarithm of the weight for that parent set size. For example, create file `weights.txt` with content

```
3
0
0
-1e100
```

This weight file specifies that we give weight 1 for parent sets of sizes 0 and 1, and an insignificantly small weight for parent sets of size 2. To sample 10 DAGs using this weight, run

```
./sampler symmetric input weights.txt 10
```

You will notice that in none of the sampled DAGs we have parent sets of size 2, as our weight for them is practically zero.

### Nonsymmetric case

In the nonsymmetric case, the weights are given in GOBNILP format:

```
<number of nodes>
for each node:
  <node name> <number of parent sets>
  for each parent set:
    <natural logarithm of weight> <number of parents> <names of parents separated by spaces>
```

The weights for parent sets that are not specified in the file are assumed to be zero. For example, create file `weights.txt` with content

```
3
A 2
0.6931471805599453 2 B C
0 0
B 2
0 0
0 1 A
C 2
0 0
0 1 A
```

For node A, the possible parent sets are {B, C} with weight 2 and the empty set with weight 1. Both B and C can either have no parents or A as a parent, both possibilities with weight 1. To sample 100 DAGs using this weight, run

```
./sampler nonsymmetric weights.txt 100
```

In the output format, the vertices are numbered in the same order as in the file, so A = 0, B = 1, C = 2.
