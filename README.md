<p align="center">
  <a href="#" title="Anton Frolov &quot;Euler's Totient&quot;"><img src="https://github.com/patnashev/prst/assets/12833063/4ec5bd6d-ebd6-443f-a3c4-7d7840fc8300"></a>
</p>

# PRST

PRST is a primality testing utility written in C++ by Pavel Atnashev. It is based on [GWnum](https://www.mersenne.org/download/) multiplication library by George Woltman.

The utility is best used for systematic searches of large prime numbers, either by public distributed projects or by private individuals. It can handle numbers of many popular forms like Proth numbers, Thabit numbers, generalized Fermat numbers, factorials, primorials and arbitrary numbers. Mersenne numbers are better handled by [GIMPS](https://www.mersenne.org/).

It is assumed that input candiates are previously sieved by a sieving utility best suited for the specific form of numbers.

The following primality tests are currently implemented:
- Fermat probabilistic test. Works for any number, although the result is not definitive.
- Proth test of k*2^n+1 numbers.
- Pocklington test of k*b^n+1 numbers, [factorial](https://en.wikipedia.org/wiki/Factorial) n!+1 and [primorial](https://en.wikipedia.org/wiki/Primorial) n#+1 numbers.
- Morrison test of k*b^n-1, n!-1, n#-1 numbers.
- LLR test of k*2^n-1 numbers, as a [special case](https://eprint.iacr.org/2023/195) of the Morrison test.

The utility can provide a high level of reliability by employing several overlapping error checks and a verification mechanism based on secure mathematical algorithms.

Donate XMR: 453MFNnC7N7NK9iK99FoNh3xZ4wY9yqt1Qk6F62iUZjg6bv6xD3oH3U2EABnUQMmco5BEK5FAhtKF18Mn43QQsysNHX9wB3,
view key to see donations: 15427a794d172323ae22fa072239e80e13a1a3f6221b43f033d460f09b167e0a.

```
Usage: PRST "<number>" <mode> <options>

Number:  K*B^N/D+C
                 K and D can occur multiple times or none.
                 K, B and D can be (<number>).
         N!+C | N!A+C | N#+C | pN#+C
                 !A is the multifactorial function.
                 pN denotes N-th prime number.
         X
                 arbitrary length number in decimal notation.
         Phi(3,[-]<number>) | Phi(6,<number>)
                 expands to X^2 ± X + 1.
         (<number>)/F
                 F can be (<number>) and occur multiple times.

Mode:    default is primality testing.
         -v
                 outputs version information.
         -info
                 outputs number profile.
         -test
                 performs built-in tests.
         -batch
                 processes several numbers.
         -order {<a> | "<number>"}
                 computes multiplicative order of a.

Options: -ini <filename>
                 reads command line options from ini file. See sample.ini.
         -log [{debug | info | warning | error}] [file <filename>]
                 level of detail of output information.
                 sets log file to write all output.
         -time [write <sec>] [progress <sec>] [coarse]
                 period of writing checkpoints and outputing progress statistics.
         -t <threads>
                 number of threads to use. For optimal performance this number
                 should depend on the amount of L3 cache.
         -spin <threads>
                 number of threads using spinwaits. Recommened values are 0 (none),
                 1 (default, only the main thread spinwaits) and the same as -t (all
                 threads spinwait, best performance but consumes 100% of CPU).
         -fft+1
         -fft [+<inc>] [safety <margin>] [generic] [info]
                 increments the size of transform used.
                 sets safety margin influencing switching to the next transform size.
                 forces generic reduction transform for debug purposes.
                 outputs only setup information, no primality test is performed.
         -cpu {SSE2 | AVX | FMA3 | AVX512F}
                 CPU instruction set to use.
         -trial
                 performs trial division by primes less than a million before the main test.
                 Slow, not a proper way to do sieving!
         -fermat [a <a>]
                 forces Fermat probabilistic test, optionally supplying starting value,
                 a = 3 by default.
         -factors [list <factor>,...] [file <filename>] [all]
                 sets list of prime factors to be used by Pocklington or Morrison tests.
                 reads the list from helper file, one factor per line.
                 forces Pocklington and Morrison tests to use all factors instead of half.
         -check [{near | always | never}] [strong [disable] [count <count>] [L <L>]]
                 enables/disables roundoff checking, by default only when close
                 to switching to the next transform size.
                 disables strong error check using Gerbicz or Gerbicz-Li algorithms.
         -proof save <count> [name <proof> <product>] [pack <name>] [keep]
                 saves proof files for verification (on tester side).
                 packs proof into a single file.
                 keeps all temporary files if asked.
         -proof build <count> [security <seed>] [roots <depth>] [name <proof> <product>]
          [pack <name>] [cert <name>] [keep]
                 builds a certificate from proof files (on authority side).
                 uses security seed to encrypt the certificate.
                 detects root of unity attack up to 2^depth degree.
                 reads proof from the pack file.
                 keeps all input files if asked.
         -proof cert {<name> | default}
                 verifies the certificate (on any side).
```
