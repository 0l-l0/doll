/+
doll - Numerical Optimization Library written in D
Copyright (C) 2019  0l-l0

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
+/

/++
Test functions for unconstrained optimization algorithms.

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.test_functions.unconstrained;

/// Himmelblau's multi-modal function. ($(LINK2 https://www.wikiwand.com/en/Himmelblau's_function,
/// Wikipedia))
real Himmelblau(in real[] x)
{
    return (x[0] ^^ 2 + x[1] - 11) ^^ 2 + (x[0] + x[1] ^^ 2 - 7) ^^ 2;
}

/// The four minimum points where Himmelblau's function takes its lowest value (0.0).
immutable HimmelblauMinPoints = [
        [3.0L, 2.0L],
        [-2.805118L, 3.131312L],
        [-3.77931L, -3.283186L],
        [3.584428L, -1.848126L]
    ];

/// Rosenbrock's banana function. ($(LINK2 https://www.wikiwand.com/en/Rosenbrock_function,
/// Wikipedia))
real Rosenbrock(real a, real b)(in real[] x)
{
    return b * (x[1] - x[0] ^^ 2) ^^ 2 + (a - x[0]) ^^ 2;
}

/// Styblinski–Tang function for multidimensional (also $(B n > 2)) tests.
real StyblinskiTang(in real[] x)
{
    import std.algorithm : map, sum;
    return x.map!"a ^^ 4 - 16 * a ^^ 2 + 5 * a".sum / 2;
}
