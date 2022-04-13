# Quantum error correction software tool
Software for finding suitable logical Z and X gates given a set of commuting Pauli generators.
### Author: Uku Kangur

This project implements the algorithm outlined by Nielsen & Chuang [[1]](#1) in Python. The project composes of four **primary** files:

<ol>
  <li>qec.py - The main program, that takes the input file, computes and writes the solution into the output file</li>
  <li>input.txt - The file that is given as input to the program</li>
  <li>outputX.txt - The file that is given as output for possible nice logical X operation variants</li>
  <li>outputZ.txt - The file that is given as output for possible nice logical Z operation variants</li>
</ol>

In addition this file has some examples of differing qubit count, that can be used to test the program.

## Starting the program

To run the program you must enter into console

```console
python qec.py input.txt
```

## input.txt format

The input Pauli group elements are given on separate lines and formatted as follows (7-qubit example):

IIIXXXX

IXXIIXX

XIXIXIX

IIIZZZZ

IZZIIZZ

ZIZIZIZ

There is no limit to the amount of lines the code can have, but the code currently supports only finding a logical X and Z for 1 logical qubit. Therefore the input set must include at least n-1 commuting generators of the stabilizer.

## Output format

The output is given in two separate files. For the logical X operation variants, the result is written into outputX.txt and for the logical Z operation variants, the result is written into outputZ.txt.

---

For any additional questions please contact the author (owner of repository).

## References
<a id="1">[1]</a> 
Nielsen, Michael A. and Chuang, Isaac L. (2010). 
Quantum Computation and Quantum Information: 10th Anniversary Edition, 425-499.
