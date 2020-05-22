# Bewley-Huggett-Aiyagari model

We compare the elapsed time for solving [the standard Aiyagari model](https://julia.quantecon.org/multi_agent_models/aiyagari.html) by each code written in Matlab, Python (with or without numba), or Julia.

We run the scripts on a machine with Intel Core i7-4790K (4.0Ghz). To measure the time, we use [hyperfine](https://github.com/sharkdp/hyperfine).

For Matlab (R2020a),
```
hyperfine 'matlab -batch "aiyagari"'
```
For Python (3.7.4 with Numba 0.45.1),
```
hyperfine 'python aiyagari.py'
```
For Julia (1.2.0),
```
hyperfine 'julia aiyagari.jl'
```

## Results

|    |  Julia  |  Python  |  Python-Numba  |  Matlab  |
| ---- | ---- | ---- | ---- | ---- |
|  Time (mean ± σ):  |  37.274 s ±  0.239 s  |    |  24.002 s ±  0.482 s  |  31.899 s ±  0.260 s  |
|  Range (min … max):  |  36.961 s … 37.758 s  |    |  23.551 s … 24.950 s  |  31.535 s … 32.215 s  |
|  Runs                |  10  |    |  10  |  10  |

## References

https://www.sas.upenn.edu/~jesusfv/Update_March_23_2018.pdf

https://web.stanford.edu/~maliars/Files/CEPR-DP13210.pdf

## Acknowledgement

I thank @yirwk and @yaoton for improving the python code.
