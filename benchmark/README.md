# SETUP

Install the necessary dependencies (Debian):

```
# apt install libbenchmark-dev libbenchmark-tools
```

# RUN A Benchmark

```
$ cmake . -B "release" -DCMAKE_BUILD_TYPE=Release
$ cd build
$ make
$ python3 /usr/share/benchmark/compare.py benchmarks BenchA BenchB
```
