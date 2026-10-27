
/**
 * @file patternMatching.c
 * @brief EPL232 Assignment 1 — DNA pattern matching with Brute Force and Karp–Rabin.
 * @author Eliana Stella Georgiou
 * @date 10/10/25
 * Reads the first line from two files (DNA text and pattern) and counts how many
 * times the pattern appears in the DNA. Supports:
 *  - "-bf"  : brute force search
 *  - "-kr"  : Karp–Rabin with base=2 and modulo Q (= INT_MAX)
 *
 * Build:
 *   gcc patternMatching.c -o patternMatching
 *
 * Run:
 *   ./patternMatching -bf dna.txt pat.txt
 *   ./patternMatching -kr dna.txt pat.txt
 */
#include <stdio.h>
#include <string.h>
#include <limits.h>

#define Q INT_MAX 
#define N 512000  /* maximum size for DNA sequence */
#define OK 0
#define ERR -1

/* FUNCTION PROTOTYPES */
/**
 * @brief Read the first line of a text file into a buffer.
 *
 * Opens @p filename for reading, reads at most @p cap-1 characters from the
 * first line into @p savedFile (leaving space for '\0'), strips a trailing
 * newline ('\n') or carriage return ('\r') if present, and closes the file.
 *
 * @param filename Path to the input text file.
 * @param savedFile Destination buffer where the line is stored
 * @param cap Capacity of @p savedFile in bytes
 * @return The length of the stored line (excluding the newline) on success;
 *         -1 on error (e.g., cannot open file, cannot read).
 */

int read_first_line(const char *filename, char *savedFile, int cap); // returns length or -1 on error

/**
 * @brief Validate that a string contains only DNA letters.
 *
 * Accepts uppercase and lowercase 'A', 'C', 'G', 'T'. Any other character
 * makes the string invalid.
 *
 * @param dna Null-terminated string to validate.
 * @return 1 if all characters are valid DNA letters; 0 otherwise.
 */

int valid_dna(const char dna[]); // 1 = valid, 0 = invalid

/**
 * @brief Parse command-line arguments: algorithm flag and two filenames.
 *
 * Expected usage:
 *   program -bf|-kr DNASequenceFile.txt patternFile.txt
 *
 * Copies the algorithm flag into @p alg and the two paths into @p dnaPath
 * and @p patPath.  minimal validation on the algorithm flag.
 *
 * @param argc Argument count from main().
 * @param argv Argument vector from main().
 * @param alg  Output buffer for "-bf" (brute force) or "-kr" (Karp–Rabin).
 * @param dnaPath Output buffer for the DNA sequence filename.
 * @param patPath Output buffer for the pattern (motif) filename.
 * @return OK (0) on success; ERR (-1) if the arguments are invalid.
 */
int parse_args(int argc, char *argv[], char alg[], char dnaPath[], char patPath[]);

/**
 * @brief Count occurrences of @p pat inside @p dna using the naive method.
 *
 * Compares @p pat against every possible starting position in @p dna and
 * counts matches (overlapping matches are allowed).
 *
 * @param dna Text string (DNA sequence).
 * @param d   Length of @p dna.
 * @param pat Pattern string to search for.
 * @param p   Length of @p pat.
 * @return Number of matches found (>= 0). Returns 0 if @p p == 0 or @p p > @p d.
 */

int brute_force(const char *dna, int d, const char *pat, int p);

/**
 * @brief Count occurrences of @p pat inside @p text using the Karp–Rabin algorithm.
 *
 * Steps:
 *  1) Precompute pow2_m_1 = 2^(m-1) % Q for the rolling update.
 *  2) Compute hash of the pattern and the first text window (length m).
 *  3) If hashes match, confirm with character-by-character comparison
 *     to avoid false positives (hash collisions).
 *  4) Slide  one position at a time:
 *       - Update the window hash by reharsh().
 *       - When hash match, confirm with compare_strings().
 *
 * Returns the total count of matches
 *
 * @param text DNA text string.
 * @param n    Length of @p text.
 * @param pat  Pattern string to search for.
 * @param m    Length of @p pat.
 * @return Number of matches found (>= 0). Returns 0 if m == 0 or m > n.
 */
long karp_rabin(const char *dna, int d, const char *pat, int p);

/**
 * @brief Compare two strings for equality up to exactly @p m characters.
 *
 * Performs a character-by-character comparison
 *
 * @param s1 Pointer to the first string segment.
 * @param s2 Pointer to the second string segment.
 * @param m  Number of characters to compare.
 * @return 1 if all @p m characters are equal; 0 otherwise.
 */

int compare_strings(const char *s1, const char *s2, int m);

/**
 * @brief Compute the initial hash of a string of length @p m (base = 2, mod Q).
 *
 * Treats the string as a base-2 polynomial in the ASCII codes:
 *   h = (s[0]*2^(m-1) + s[1]*2^(m-2) + ... + s[m-1]*2^0) mod Q
 *
 * Implemented iteratively as:
 *   h = 0;
 *   for i in [0..m-1]: h = (h*2 + (unsigned char)s[i]) % Q;
 *
 * Using modulo at each step keeps h bounded in [0, Q-1].
 *
 * @param s Input string (pattern or the first window of text).
 * @param m Number of characters to include in the hash.
 * @return The computed hash value as a signed 64-bit integer (long long).
 */
long long hash(const char *s, int m);

/**
 * @brief Rolling-hash update when sliding the window by one to the right.
 *
 * Given the previous window hash @p h for text[i..i+m-1], update the hash for
 * the next text[i+1..i+m] using:
 *   h' = ((h - left*2^(m-1)) * 2 + right) mod Q
 *
 * Where:
 *   - @p left  is the character that leaves the window (text[i])
 *   - @p right is the character that enters the window (text[i+m])
 *   - @p pow2 is 2^(m-1) mod Q (precomputed)
 *
 * A correction is applied if the intermediate subtraction becomes negative,
 * by adding Q before continuing.
 *
 * @param h Current hash of the previous window.
 * @param left Unsigned value of the character leaving the square.
 * @param right Unsigned value of the character entering the square.
 * @param pow2 Precomputed (2^(m-1)) % Q.
 * @return The new  hash value modulo Q.
 */
long long reharsh(long long h, unsigned left, unsigned right, long long pow2);

/**
 * @brief Entry point of the program.
 *
 * Handles command-line arguments, reads input files, validates DNA sequences,
 * and calls the selected search algorithm (Brute Force or Karp–Rabin).
 *
 * Usage example:
 *   ./patternMatching -bf dnaSequence.txt patSequence.txt
 *   ./patternMatching -kr dnaSequence.txt patSequence.txt
 *
 * @param argc Number of command-line arguments (should be 4).
 * @param argv Array of argument strings:
 *             - argv[0]: program name
 *             - argv[1]: algorithm flag ("-bf" or "-kr")
 *             - argv[2]: DNA sequence file
 *             - argv[3]: pattern file
 * @return 0 on successful execution; non-zero on error.
 */
int main(int argc, char *argv[])
{
    char alg[4];                     
    char dnaPath[260], patPath[260]; /* filenames for DNA and pattern */
    static char dnaSeq[N + 1];       /* static arrays to avoid stack overflow   */
    static char patSeq[N + 1];
    int dnaLen = 0, patLen = 0; /* lengths of sequences */

    /* parse command line arguments */
    if (parse_args(argc, argv, alg, dnaPath, patPath) == ERR)
    {
        printf("Error you missed an argument! :( \n");
        return 1;
    }

    /* read DNA sequence (returns length or -1) */
    dnaLen = read_first_line(dnaPath, dnaSeq, N + 1);
    if (dnaLen < 0)
        return 1;

    /* read pattern sequence */
    patLen = read_first_line(patPath, patSeq, N + 1);
    if (patLen < 0)
        return 1;

    /* check sequence lengths against maximum N */
    if (dnaLen >= N)
    {
        printf("DNA sequence too large\n");
        return 1;
    }
    if (patLen >= N)
    {
        printf("Pattern sequence too large\n");
        return 1;
    }

    /* validate sequences: must contain only A, C, G, T (uppercase or lowercase) */
    if (!valid_dna(dnaSeq))
    {
        printf("Invalid characters in DNA sequence (allowed: A,C,G,T)\n");
        return 1;
    }
    if (!valid_dna(patSeq))
    {
        printf("Invalid characters in pattern (allowed: A,C,G,T)\n");
        return 1;
    }

    int count = 0;
    if (strcmp(alg, "-bf") == 0)
    {
        /* run brute force algorithm */
        count = brute_force(dnaSeq, dnaLen, patSeq, patLen);
    }
    else
    {
        /* run Karp-Rabin algorithm */
        count = karp_rabin(dnaSeq, dnaLen, patSeq, patLen);
    }

    /* print results */
    printf("The pattern was found: %d times\n :)", count);

    return 0;
}

int read_first_line(const char *filename, char *savedFile, int cap)
{
    FILE *f = fopen(filename, "r");
    if (!f)
    {
        printf("Cannot open file: %s\n", filename);
        return -1;
    }

    if (!fgets(savedFile, cap, f)) // read line into buffer
    {
        printf("Empty file or read error: %s\n", filename);
        fclose(f);
        return -1;
    }
    fclose(f);

    int len = (int)strlen(savedFile);
    /* remove newline or carriage return at end if present */
    while (len > 0 && (savedFile[len - 1] == '\n' || savedFile[len - 1] == '\r'))
    {
        savedFile[len - 1] = '\0';
        len--;
    }
    return len;
}

int valid_dna(const char dna[])
{
    int i = 0;
    while (dna[i] != '\0')
    {
        char c = dna[i];
        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
              c == 'a' || c == 'c' || c == 'g' || c == 't'))
        {
            return 0; // invalid character found
        }
        i++;
    }
    return 1;
}

int parse_args(int argc, char *argv[], char alg[], char dnaPath[], char patPath[])
{
    if (argc != 4)
        return ERR;

    strncpy(alg, argv[1], 3); // copy algorithm flag
    alg[3] = '\0';
    strncpy(dnaPath, argv[2], 259); // copy DNA filename
    dnaPath[259] = '\0';
    strncpy(patPath, argv[3], 259); // copy pattern filename
    patPath[259] = '\0';

    if (strcmp(alg, "-bf") != 0 && strcmp(alg, "-kr") != 0)
        return ERR;

    return OK;
}

int brute_force(const char *dna, int d, const char *pat, int p)
{
    if (p > d || p == 0)
        return 0; /* pattern cannot fit */

    int found = 0;

    /* check each possible start position */
    for (int i = 0; i <= d - p; ++i)
    {
        int j;
        for (j = 0; j < p; ++j)
        {
            if (pat[j] != dna[i + j])
            {
                break; /* exit loop */
            }
        }
        if (j == p)
        { /* all characters matched */
            ++found;
        }
    }

    return found;
}

int compare_strings(const char *s1, const char *s2, int m)
{
    for (int i = 0; i < m; i++)
    {
        if (s1[i] != s2[i])
            return 0;
    }
    return 1;
}

long long hash(const char *s, int m)
{
    long long h = 0;
    for (int i = 0; i < m; i++)
    {
        h = (h * 2 + (unsigned char)s[i]) % Q;
    }
    return h;
}

long long reharsh(long long h, unsigned left, unsigned right, long long pow2)
{
    h = (h - (long long)left * pow2) % Q; // remove contribution of left char
    if (h < 0)
        h += Q;              // fix negative remainder
    h = (h * 2 + right) % Q; 
    return h;
}

long karp_rabin(const char *text, int n, const char *pat, int m)
{
    if (m == 0 || m > n)
        return 0;

    /* compute pow2_m_1 = 2^(m-1) mod Q */
    long long pow2 = 1LL;
    for (int i = 0; i < m - 1; i++)
    {
        pow2 = (pow2 * 2) % Q;
    }

    long long hP = hash(pat, m);  /* hash of pattern */
    long long hS = hash(text, m); /* hash of first sequence */

    long count = 0;

    /*  first window */
    if (hS == hP && compare_strings(text, pat, m))
    {
        ++count;
    }

    /* slide window */
    for (int i = 1; i + m <= n; i++)
    {
        hS = reharsh(hS,
                     (unsigned char)text[i - 1],
                     (unsigned char)text[i + m - 1],
                     pow2);

        if (hS == hP && compare_strings(text + i, pat, m))
        {
            ++count;
        }
    }
    return count;
}
