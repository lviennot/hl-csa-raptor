# hl-csa-raptor
Hub labeling based version of CSA and RAPTOR supporting unrestricted walking

An implementation in C++ of HLCSA and HLRaptor which are variants CSA and RAPTOR supporting unrestricted walking through the use of a hub labeling of the footpath graph.

Associated data can be found at this [graph repository](https://files.inria.fr/gang/graphs/public_transport/)

Supported queries:
 * Earliest arrival time.
 * Multi-criteria (arrival time, number of transfers, walking time), RAPTOR only.
 * Profile (a.k.a. range).

## Compile:
```
make
```

## Try it (with download of data):
```
make test
```

Author : Laurent Viennot, Inria 2019

Acknowledgement to Duc-Minh Phan for implementing a preliminary version of the CSA part.  
