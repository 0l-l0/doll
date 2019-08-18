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
Often used test functions to test the convergence, robustness and precision of
an optimization algorithm. Two collections of these functions are available in
this package:
$(UL
    $(LI Constrained)
    $(LI Unconstrained)
)
An exhaustive collection of these kind of functions can be found
$(LINK2 https://www.wikiwand.com/en/Test_functions_for_optimization, here).

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.test_functions;
