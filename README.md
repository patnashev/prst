# PRST
PRST is a primality testing utility written in C++.

```
Usage: PRST {"K*B^N+C" | <file>} <options>
Options: [-t <threads>] [-spin <threads>] [-fft+1] [-log {debug | info | warning | error}] [-time [write <sec>] [progress <sec>]]
         -fermat [a <a>]
         -proof {save <count> | build <count> [security <seed>] | cert {<name> | default}} [name <proof> <product> [{<cert> | default}]]
         -check [{near | always | never}] [strong [count <count>] [L <L>]]
```
