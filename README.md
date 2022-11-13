# PRST
PRST is a primality testing utility written in C++.

```
Usage: PRST {"K*B^N+C" | "N!+C" | "N#+C" | "N"} <options>
Options: -log {debug | info | warning | error}
         -t <threads>
         -spin <threads>
         -time [write <sec>] [progress <sec>]
         -fft+1
         -fft [+<inc>] [safety <margin>]
         -cpu {SSE2 | AVX | FMA3 | AVX512F}
         -fermat [a <a>]
         -proof {save <count> | build <count> [security <seed>] [roots <depth>] | cert {<name> | default}} [name <proof> <product> [{<cert> | default}]]
         -check [{near | always| never}] [strong [count <count>] [L <L>]]
```
