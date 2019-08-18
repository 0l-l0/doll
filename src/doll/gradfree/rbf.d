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
Radial basis function optimizer. Useful resources:
$(UL
    $(LI doll/literature/RBFOpt_Library_Paper.pdf)
    $(LI doll/literature/RBF_Matrix_Factorization_Techniques.pdf)
    $(LI https://github.com/coin-or/rbfopt)
    $(LI https://developer.ibm.com/open/openprojects/rbfopt/)
)

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.gradfree.rbf;

version (WithLapack):

import doll.optimizer_skeletons : RealOptimizer;
import doll.utils : funRnToR1_t, Bounds;

// ########## COMMON RBF FUNCTIONS ##########

/// Linear function (identity).
real linear(real x)
{
    return x;
}

/// Cubic.
real cubic(real x)
{
    return x ^^ 3;
}

/// Multiquadric. `gamma` is a shape parameter usually set to 1.
real multiquadric(real gamma)(real x)
{
    return sqrt(x ^^ 2 + gamma ^^ 2);
}

/// Thin plate spline function.
real thinPlateSpline(real x)
{
    import std.math : log;
    return x ^^ 2 * log(x);
}

// ########## OPTIMIZER CLASS ##########

/++
$(B $(RED NOT IMPLEMENTED YET!))
+/
class RBF : RealOptimizer
{
    real[][] pointSet; // Why should it be public? It shouldn't.
    private
    {
        size_t _dim; // Replace it with a local 'dim' in 'minimize()' (Think it over!)
    }

    // CONSTRUCTORS

    // Initialization without bounds is not allowed (Or should it be?)
    //mixin(optimizerBaseConstructors);

    ///
    this(funRnToR1_t objFun_, real[][2] bounds_)
    {
        this(objFun_, Bounds(bounds_[0], bounds_[1]));
    }

    ///
    this(funRnToR1_t objFun_, Bounds bounds_)
    {
        super(objFun_);
        bounds = bounds_;
        _dim = bounds_.dim;
    }

    // MEMBER FUNCTIONS

    // It is very sketchy and does not minimize anything at all!
    override void minimize()
    {
        import doll.utils : euclideanNorm;
        alias RBFFun = thinPlateSpline;
        immutable minPolyDim = 1; // For thin plate spline
        immutable diagValue = 0.0; // RBFFun(0);
        immutable initialPointNum = _dim + 1;

        initPointSet;

        real[][] interpMatrix;
        interpMatrix.length = initialPointNum + minPolyDim * _dim + 1; // Usually: 2 * _dim + 2
        foreach (ref row; interpMatrix)
            row.length = interpMatrix.length;
        real[] diff;
        diff.length = _dim;
        foreach (pointInd1, point1; pointSet)
        {
            interpMatrix[pointInd1][pointInd1] = diagValue;
            foreach (pointInd2, point2; pointSet[pointInd1 + 1 .. $])
            {
                diff[] = point1[] - point2[];
                auto adjustedInd2 = pointInd2 + pointInd1 + 1;
                interpMatrix[pointInd1][adjustedInd2] =
                interpMatrix[adjustedInd2][pointInd1] = RBFFun(diff.euclideanNorm);
            }
            interpMatrix[pointInd1][initialPointNum .. $] = point1[] ~ 1.0; // Is it efficient enough? (~)
            interpMatrix[pointInd1 + initialPointNum][initialPointNum .. $] = .0;
            foreach (rowInd, ref row; interpMatrix[initialPointNum .. $ - 1])
                row[pointInd1] = point1[rowInd];
        }
        interpMatrix[$ - 1][0 .. initialPointNum][] = 1.0;
    }

    private void initPointSet(size_t numPoints = 0)
    {
        if (1 == _dim)
        {
            pointSet.length = 2;
            pointSet[0] = bounds.lower.dup;
            pointSet[1] = bounds.upper.dup;
        }
        else
        {
            // If no 'numPoints' value is given (i.e. set to zero), we create
            // '_dim + 1' points
            if (!numPoints)
                numPoints = _dim + 1;

            // TODO: Reduce the local imports if possible (they are better than
            // module level "global" imports but if it is feasible, pull in the
            // whole module at top level)
            import doll.utils : eigenValues, latinHypercubeDesignMaxiMin;
            import std.algorithm : minElement;

            immutable maxTrials = 55;
            immutable dependencyThreshold = 1.0e-6;

            real[][] LatinHypercube;

            foreach (_; 0 .. maxTrials)
            {
                LatinHypercube = latinHypercubeDesignMaxiMin(bounds, numPoints); // Does it work? (assigment)
                auto eigenValsOfLH = eigenValues(cast(double[][]) LatinHypercube);

                if (dependencyThreshold < eigenValsOfLH.minElement)
                    break;
            }

            // TODO: Check whether it has been succeeded to find a feasible LH
            // in 'maxTrials' iterations. If not, handle it. (E.g. exception or
            // some not feasible default LH etc.)

            pointSet.length = numPoints;
            foreach (pointInd, ref point; pointSet)
                point = LatinHypercube[pointInd]; //.dup ?
        }
    }
}

/+
unittest
{}
+/
