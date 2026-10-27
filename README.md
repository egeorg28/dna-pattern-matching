# DNA Pattern Matching

A command-line program written in C that counts occurrences of a DNA pattern within a sequence using Brute Force and Rabin–Karp algorithms.

Developed for the EPL232 course at the University of Cyprus.

## Features

- Brute Force pattern matching.
- Rabin–Karp matching using a rolling hash.
- Character-by-character verification of hash matches to handle collisions.
- Detection of overlapping matches.
- Input validation for DNA characters: A, C, G and T, in uppercase or lowercase.
- File input and command-line algorithm selection.

## Requirements

- GCC or another C compiler supporting C99 or later.

## Build

```bash
gcc -std=c11 -Wall -Wextra -Wpedantic patternMatching.c -o patternMatching
```

## Usage

```bash
./patternMatching <algorithm> <dna-file> <pattern-file>
```

Available algorithms:

- `-bf`: Brute Force
- `-kr`: Rabin–Karp

## Example

Run either algorithm with the included sample files:

```bash
./patternMatching -bf dnaSequence.txt patSequence.txt
./patternMatching -kr dnaSequence.txt patSequence.txt
```

Both algorithms find **3 occurrences** of the pattern `gaa`.

## Implementation

Brute Force compares the pattern at every possible starting position.

Rabin–Karp uses a rolling hash with base 2 and modulus `INT_MAX`.
When hashes match, the program compares the characters to confirm the match.

Hash calculations use `long long` intermediates to accommodate multiplication before applying the modulus.

## Complexity

For a sequence of length n and a pattern of length m:

- **Brute Force:** O(nm) worst-case time.
- **Rabin–Karp:** O(n + m) hashing work, plus O(m) for each hash-match verification; O(nm) worst-case time.

Actual performance depends on the input, pattern length and hash collisions.

## Input Behaviour and Limitations

- Only the first line of each input file is processed.
- Matching is case-sensitive: `ACG` and `acg` are treated as different strings.
- Invalid DNA characters cause an error.
- An empty pattern or a pattern longer than the sequence returns zero matches.
- Input buffers have a fixed capacity. The current implementation rejects sequences of 512,000 characters or more.

## Author

Eliana Stella Georgiou
