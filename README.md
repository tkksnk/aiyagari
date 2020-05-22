# Bewley-Huggett-Aiyagari model

We compare the elapsed time for each code written in Matlab, Python (with or without numba), or Julia.

We run the scripts on a machine with Intel Core i7-4790K (4.0Ghz). To measure the time, we use hyperfine.

For Matlab,
```
hyperfine 'matlab -batch "aiyagari"'
```
For Python,
```
hyperfine 'python aiyagari.py'
```
For Julia,
```
hyperfine 'julia aiyagari.jl'
```

## Results

|  Language  |  Matlab  |  Python  |  Python-Numba  |  Julia  |
| ---- | ---- | ---- | ---- | ---- |
|  Time (mean ± σ):  |  37.274 s ±  0.239 s  |    |  TD  |  TD  |
|  Range (min … max):  |  36.961 s … 37.758 s  |    |  TD  |  TD  |
