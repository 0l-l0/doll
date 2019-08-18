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
First of all, this module contains the optimizer class implementin the SPSA
(Simultanious Perturbation Stochastic Approximation) algorithm. This method is a
kind of stochastic optimization algorithm which operates like some of those
methods using finite-differences but unlike them it computes only two
measurement of the objective funtction in each step. (Usually, the others work
with $(B 2n) measurements.)
For more details see:
$(UL
    $(LI [1] http://www.jhuapl.edu/SPSA/)
    $(LI [2] doll/literature/Spall_Implementation_of_SPSA_1998.pdf)
    $(LI [3] doll/literature/Spall_Stochastic_Optimization_2012)
    $(LI [4] doll/literature/Adaptive_Init_Step_Size_Selection_for_SPSA_2016)
)
Furthermore, we define here a parameter-tuple to facilitate the initialization
of the actual optimizer instance.

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.gradfree.spsa;

import doll.optimizer_skeletons : optimizerBaseConstructors, RealOptimizer;
import doll.utils : funRnToR1_t;

import std.typecons : Tuple;

/// SPSA parameter type to wrap up the typical parameters of the algorithm.
alias SPSAParameters = Tuple!(
    real, "alpha",
    real, "gamma",
    real, "A",
    real, "a",
    real, "c"
);

/++
Class for the SPSA (Simultanious Perturbation Stochastic Approximation) algorithm.
+/
class SPSA : RealOptimizer
{
    /// Default values are chosen after Spall's recommendation [2] but further
    /// improvement is achievable with adaptive step-size implementations (e.g. [4]).
    auto param = SPSAParameters(
        .602, // alpha (based on Spall's recommendation [2])
        .101, // gamma (based on Spall's recommendation [2])
        100,  // A (chosen to be floor(maxIter / 10))
        .06,  // a (tested on Himmelblau's function)
        .76   // c (tested on Himmelblau's function)
    );
    
    /// $(B n)-dimensional starting point for the algorithm.
    real[] thetaInit;
    
    // CONSTRUCTORS

    mixin(optimizerBaseConstructors);
    
    /// Constructor with the given objective function and an initial theta value
    /// where SPSA should start from. If we set bounds for the search process,
    /// it is recommended to choose `thetaInit` from the hypercube defined by
    /// the two endpoints.
    this(funRnToR1_t objFun_, real[] thetaInit_)
    {
        super(objFun_);
        thetaInit = thetaInit_.dup;
    }

    /// We can construct an SPSA object giving the objective function
    /// with a predefined (e.g. read from file) parameter tuple.
    this(funRnToR1_t objFun_, SPSAParameters param_)
    {
        super(objFun_);
        param = param_;
    }

    /// Constructor for the case when we set all SPSA specific initial values.
    this(funRnToR1_t objFun_, real[] thetaInit_, SPSAParameters param_)
    {
        this(objFun_, thetaInit_);
        param = param_;
    }

    // MEMBER FUNCTIONS

    override void minimize()
    {
        immutable dim = thetaInit.length;

        // Without 'objFun' or given starting point it does nothing
        if (null == objFun || size_t.init == dim)
            return;

        // INIT SECTION

        import std.algorithm : max, min;
        import std.random : choice;

        real[] theta, thetaPlus, thetaMinus, delta, grad;
        theta.length      =
        thetaPlus.length  =
        thetaMinus.length =
        delta.length      =
        grad.length       = dim;

        theta = thetaInit.dup;

        real ak, ck;
        
        // MAIN LOOP

        foreach (k; 1 .. maxIter)
        {
            // Update of the gain sequences
            ak = param.a / (k + param.A) ^^ param.alpha;
            ck = param.c / k ^^ param.gamma;

            // Generation of the Bernoulli vector and using it to
            // construct the perturbations of the previous theta value
            foreach (ref coord; delta)
                coord = choice([-1, 1]);
            thetaPlus[] = theta[] + ck * delta[];
            thetaMinus[] = theta[] - ck * delta[];

            // Computation of the k-th gradient approximation
            grad[] = (objFun(thetaPlus) - objFun(thetaMinus)) /
                (2 * ck * delta[]);
            
            // Move towards the optimum
            theta[] -= ak * grad[];

            // If bounded, keep 'theta' between the bounds
            if (!bounds.empty)
                foreach (ind, ref th; theta)
                    th = th.max(bounds.lower[ind])
                           .min(bounds.upper[ind]);
        }

        // SET MINIMUM VALUE AND POINT
        //minimumPoint.length = dim;
        minimumPoint = theta;
        minimumValue = objFun(theta);
    }
}

/// As we can see in the following example, the precise choice of the parameters
/// is very important due to the statistical sensitivity of SPSA.
unittest
{
    import doll.test_functions.unconstrained;
    import std.algorithm : any, map, minIndex;
    import std.math : approxEqual;
    import std.random : uniform;

    // These parameters seem really good for Himmelblau's function
    immutable param = SPSAParameters(.5, .5, 90, .06, .76);
	auto SPSAOpt = new SPSA(&Himmelblau, param);

    // In the four iteration steps we will visit all quadrants of the
    // two dimensional plane (the bounds specify the current quadrant)
    real[2][2] bounds;

    foreach (k; 0 .. 4)
    {
        // Some tricky quadrant calculation
        bounds[0] = -4.0L * [k >> 1, k & 1];
        bounds[1] = 4.0L * [(3 - k) >> 1, (k + 1) & 1];

        SPSAOpt.bounds = bounds;

        // We execute ten independent minimization processes one after another
        // in the same quadrant but with different randomly-chosen starting points
        real[][] minPoints;
        real[] minVals;

        foreach (_; 0 .. 10)
        {
            SPSAOpt.thetaInit = [
                uniform(bounds[0][0], bounds[1][0]),
                uniform(bounds[0][1], bounds[1][1])
            ];

            SPSAOpt.minimize;
            minPoints ~= SPSAOpt.minPoint;
            minVals ~= SPSAOpt.minVal;
        }

        // Selection of the point with the lowest objective function value
        auto bestMinPoint = minPoints[minVals.minIndex];

        // We hope it is close enough to one of the known minimum points of
        // Himmelblau's function
        assert(HimmelblauMinPoints.map!(a => approxEqual(a, bestMinPoint)).any);
    }
}
