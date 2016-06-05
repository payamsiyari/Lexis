# Lexis
Lexis: An Optimization Framework for Discovering the Hierarchical Structure of Sequential Data

This is the implementation of the Lexis framework presented in the paper: http://arxiv.org/abs/1602.05561

## Requirements
* The code uses a slightly modified (in terms of I/O) version of the suffix array library, proposed in [1].
* The code uses NetworkX package [2] for graph algorithms.

## Installation
* Make sure to compile the code in repeats1 directory by simply running ```make```.
* Make sure to install NetworkX by running ```pip install networkx```.

## Usage
```python
./python Lexis.py [-t (c | i | s) | -p (i) | -q | -r (r | mr | lmr | smr) | -f (c | e) | -m | -l] <filename>
    [-t]: choosing between character sequence, integer sequence or space-separated sequence
        c - character sequence
        i - integer sequence
        s - space-separated sequence
    [-p]: specifies DAG printing option (for debugging purposes)
        i - prints the DAG in integer sequence format
    [-q]: disables logging
    [-r]: repeat type (for normal repeat replacements)
        r - repeat
        mr - maximal repeat (default)
        lmr - largest-maximal repeat
        smr - super-maximal repeat
    [-f]: cost function to be optimized
        c - concatenation cost
        e - edge cost (default)
    [-m]: consider each line of the input file as a separate target string
    [-l]: load a DAG file (will override -r -t -m options)
```

## References
* [1] M. Gallé. ''Searching for compact hierarchical structures in DNA by means of the Smallest Grammar Problem'', PhD thesis, Université Rennes 1, 2011.
* [2] https://networkx.github.io