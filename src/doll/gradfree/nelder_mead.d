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
The downhill simplex or amoeba method is one of the simplest, yet very effective
algorithms. This method uses only $(B n+1) initial measurements of the objective
function, and on average it updates one of these points in every iteration (in
the worst case it updates $(B n) of them). This method is often named after its
inventors and we will also refer to it as
$(LINK2 https://www.wikiwand.com/en/Nelder%E2%80%93Mead_method, Nelder–Mead).

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.gradfree.nelder_mead;

import doll.optimizer_skeletons : RealOptimizer, optimizerBaseConstructors;
import doll.utils : funRnToR1_t, Simplex;

/++
The amoeba method implementation.
+/
class NelderMead : RealOptimizer
{
    /// Simplex object to store the actual "sampling points" in each step.
    Simplex simplex;

    /// A counter to query the actually needed iteration steps after running
    /// `minimize()`.
    uint neededIter;

    // CONSTRUCTORS
    
    mixin(optimizerBaseConstructors);

    /// Constructor with an objective function and an initial simplex which
    /// might be chosen with respect to some heuristic or in completely random manner.
    this(funRnToR1_t objFun_, Simplex simplex_)
    {
        super(objFun_);
        simplex = simplex_;
    }

    // MEMBER FUNCTIONS

    // TODO: Implement bound checking for the case when we set the 'bounds' member
    override void minimize()
    {
        immutable dim = simplex.dim;

        // Without 'objFun' or given initial simplex it does nothing
        if (null == objFun || size_t.init == dim)
            return;
        
        // INIT SECTION

        import doll.utils : euclideanNorm;

        import std.algorithm : each, map, minIndex, sum, topNIndex;
        import std.conv : to;
        import std.math : cmp, sqrt;
        import std.typecons : Yes;

        immutable reflDistances = [2, 1.7, 1.4, 1, .7, -.5, .4, -.3];
        real[reflDistances.length] candidateFunValues;
        real[][reflDistances.length] newPointCandidates;

        size_t[2] highestIndices;
        size_t lowestIndex;

        real[] center, reflContDirection, tmp;
        reflContDirection.length =
                      tmp.length =
                   center.length = dim;

        // Objective function values of the initial simplex points
        //auto funVals = simplex.points.map!objFun.array; // with comp.-time 'objFun'
        real[] funVals;
        funVals.length = dim + 1;
        foreach (ind, ref fx; funVals)
            fx = objFun(simplex.points[ind]);

        // MAIN LOOP

        neededIter = maxIter;
        foreach (iter; 0 .. maxIter)
        {
            // Select the two points with the highest 'objFun' values
            topNIndex!((a, b) => cmp(a, b) > 0)
                (funVals, highestIndices[], Yes.sortOutput);

            // Reflection of the point with the highest 'objFun' value through
            // the centroid of the other simplex points (all points but itself)
            center = simplex.centroid(highestIndices[0]);
            reflContDirection[] = center[] - simplex.points[highestIndices[0]][];

            foreach (ind, ref ptCand; newPointCandidates)
            {
                tmp[] = center[] + reflDistances[ind] * reflContDirection[];
                ptCand = tmp.dup;
            }
            // If we knew the 'objFun' at compile-time, we could
            // simplify the following loop this way:
            //candidateFunValues = newPointCandidates.map!objFun.array;
            foreach (ind, ref fx; candidateFunValues)
                fx = objFun(newPointCandidates[ind]);
            auto candMinInd = candidateFunValues[].minIndex;

            // After we choose the candidate with minimal 'objFun' value
            // we check its feasibility (i.e. has it lower 'objFun' value than
            // the point with the second highest value)
            // If not, we contract the simplex towards the 'lowestPoint' (the
            // one with the lowest 'objFun' value)
            if (candidateFunValues[candMinInd] >= funVals[highestIndices[1]])
            {
                lowestIndex = funVals.minIndex;
                auto lowestPoint = simplex.points[lowestIndex];

                foreach (ind, ref pt; simplex.points)
                    if (lowestIndex != ind)
                    {
                        tmp[] = lowestPoint[] - pt[];
                        pt[] += .5 * tmp[];
                        funVals[ind] = objFun(pt);

                        reflContDirection[] += tmp[];
                    }
                reflContDirection[] /= dim;
            }
            // If the candidate is feasible, we will replace the simplex point
            // with the previously measured highest 'objFun' value with the new point
            else
            {
                simplex.points[highestIndices[0]] = newPointCandidates[candMinInd];
                funVals[highestIndices[0]] = candidateFunValues[candMinInd];
            }

            // TODO: Introduce some other (better) optimality checks
            if (reflContDirection.euclideanNorm
                .cmp(real.min_normal.sqrt * 4) < 0) // Less than the sqrt(mach.prec.) * 4
            {
                neededIter = iter;
                break;
            }
        }
        // Storage of the minimum point and value
        minimumPoint.length = dim;
        minimumPoint[] = simplex.pointSum[] / (dim + 1);
        minimumValue = objFun(minimumPoint);
    }
}

///
unittest
{
    import doll.test_functions.unconstrained;
    import std.algorithm : any, map;
    import std.math : approxEqual, cmp;

    // Testing with a randomly initialized two-dimensional simplex and
    // Rosenbrock's test function (its minimum is at: (11, 11^2))
    auto simplex = Simplex(2);
    auto NMOpt = new NelderMead(&Rosenbrock!(11,100), simplex);

    NMOpt.minimize;
    assert(approxEqual(NMOpt.simplex.points, [[11, 121], [11, 121], [11, 121]]));
    assert(approxEqual(NMOpt.minPoint, [11.0L, 121.0L]));
    assert(cmp(NMOpt.minVal, 1e-31) < 0); // Hopefully it is below 1e-31

    // Reinitializing the simplex object and changing the objective function
    simplex.randomInit(2);
	NMOpt.simplex = simplex;
	NMOpt.objFun = &Himmelblau;

    // Checking all the possible minimum points of Himmelblau's test function
	NMOpt.minimize;
    assert(HimmelblauMinPoints.map!(a => approxEqual(a, NMOpt.minPoint)).any);

    // Without any given objective function or simplex it does not minimize anything
    auto emptyNelderMead = new NelderMead();
    emptyNelderMead.minimize;
    assert(emptyNelderMead.minPoint == []);
    assert(!cmp(emptyNelderMead.minVal, real.init));
    assert(emptyNelderMead.neededIter == 0);
}
