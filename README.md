# NWPA_Cpp
Implementation of my favorite Needleman-Wunsch pairwise alignment algorithm. Based on my NWPA package for R language.

The `align` function is the main function of this snippet. It accepts a couple of nucleotide sequences (gapped or with ambiguous base pairs) and returns a tuple with two aligned sequences. For example, the code 
```
std::string al1, al2;
std::tie(al1, al2) = align("GCATGCU", "GATTACA", 1, -1, -1);
std::cout << al1 << std::endl;
std::cout << al2 << std::endl;
```
Will output 
```
GCA-TGCU
G-ATTACA
```
