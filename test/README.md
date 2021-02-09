# SETUP

Install the necessary dependencies:

```
# apt install googletest googletest-tools cmake
# sudo apt-get install libgmp3-dev
# cd /usr/src/googletest
# cmake .
# sudo cmake --build . --target install
```

# RUN TESTS

```
$ mkdir tests
$ cd tests
$ cmake ..
$ make
$ make tests
```
