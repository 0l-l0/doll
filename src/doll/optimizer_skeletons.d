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
Abstract base classes encapsulating the common objects and operations for any
future optimizers. This module also contains some commonly used code snippets
in form of template mixins, string mixins and even
$(LINK2 https://dlang.org/blog/2017/08/28/open-methods-from-c-to-d/, open methods).

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.optimizer_skeletons;

import doll.utils : funRnToR1_t, Bounds;

/// Simple mixin code for base constructors.
enum optimizerBaseConstructors =
q{
    /// Empty constructor.
    this() { super(); }
    /// Construct the optimizer only with the objective function to be minimized.
    this(funRnToR1_t objFun_) { super(objFun_); }
};

/+ /// An extendable version of the base constructor mixin string
template optimizerBaseConstructors(string commonCode)
{
    const char[] optimizerBaseConstructors =
    "
        /// Empty constructor.
        this() { super();" ~ commonCode ~ "}
        /// Construct the optimizer only with the objective function to be minimized.
        this(funRnToR1_t objFun_) { super(objFun_);" ~ commonCode ~ "}
    ";
}
+/

/++
Base class for all of the optimizers implementing a minimization algorithm on
some $(B R^n -> R) funtion.
+/
abstract class RealOptimizer
{
    /// Objective function. This delegate stores the function to
    /// be minimized. Only functions with `real function(in real[])`
    /// ($(B R^n -> R)) are allowed.
    funRnToR1_t objFun;

    /// Maximum number of iterations of the main loop in the actual algorithm.
    /// Be careful, it may differ from the number of real function evaluations!
    /// Default value is 1000.
    uint maxIter = 1000;

    /// A dimension-checked `Bounds` variable containing two $(B n)-dimensional
    /// point values which define a hypercube to perform the optimization
    /// process in it.
    Bounds bounds;

    protected
    {
        real minimumValue;
        real[] minimumPoint;
    }
    
    /// Empty constructor.
    this() {}
    
    /// Construct the optimizer only with the objective function to be minimized.
    this(funRnToR1_t objFun_)
    {
        objFun = objFun_;
    }

    /// Getter property for querying the minumum value computed the by the last
    /// run of `minimize()`.
    final @property real minVal() { return minimumValue; }

    /// Getter property for querying the minumum point computed the by the last
    /// run of `minimize()`.
    final @property real[] minPoint() { return minimumPoint.dup; }

    /// Performs the actual optimization process. This function wraps the
    /// different algorithms implemented by descendants of `RealOptimizer`.
    abstract void minimize();
}
