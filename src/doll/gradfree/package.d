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
These algorithms are very diverse in the sense of the mathematical concepts they
are built on but they share a common trait, namely the independency from the
gradient of the objective function.
$(BR)Implemented optimizers:
$(UL
    $(LI Nelder–Mead)
    $(LI RBF)
    $(LI SPSA)
)

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.gradfree;

public
{
    import doll.gradfree.nelder_mead;
    import doll.gradfree.rbf;
    import doll.gradfree.spsa;
}
