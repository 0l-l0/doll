# doll
[![Build Status](https://travis-ci.com/0l-l0/doll.svg?branch=master)](https://travis-ci.com/0l-l0/doll)
[![codecov](https://codecov.io/gh/0l-l0/doll/branch/master/graph/badge.svg)](https://codecov.io/gh/0l-l0/doll)

**D** **O**ptimization **L**ibrary for **L**earning - A _doll_ you can play with,
dress it up and imitate some real-life situations, as if it was alive.

doll is a collection of several gradient-free and gradient-dependent optimization
algorithms. The primary goal of this library was to test the usability of the
[D language](https://dlang.org/) for such mathematical purposes.

## Documentation
For API reference you can generate the library documentation from the source
code using one of the D documentation generators (ddoc, ddox, adrdoc etc.). For
example run:

```sh
dub build -b ddox
```

## Testing
To test the implemented optimization functions you can write your own test
project or use the example codes from the documentation. These are also unit
test sections therefore you can easily test all of them with the following
command.

```sh
dub test
```

Or if you are using ldc:

```sh
dub test --compiler=ldc2
```

## TODO list

* [ ] Implement [ADAM](https://arxiv.org/pdf/1412.6980.pdf) optimizer (gradient dependent)
* [ ] Optional parallelism in certain calculations.
* [ ] OTHER: Use D's recurrence and sequence methods (and ranges in general) to
express mathematical formulas

## License
The library is distributed under the terms and conditions of GNU LGPL v3.0. 
Copyright &copy; 2019, 0l-l0

