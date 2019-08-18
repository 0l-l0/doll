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
Utilities, type aliases.

Authors: 0l-l0
License: LGPL-3.0
Copyright: 0l-l0
+/
module doll.utils;

/+ For CTFE purposes
+  alias Point(uint dim) = double[dim];
+/

// ########## ALIASES ##########

/++
Type alias for general objective functions.
[$(I $(D_KEYWORD in) qualifier may be too strict in some cases of more
complicated functions because of its `scope` qualifier wrapped up by
it (which currently is experimental but in the future it may impact on the
behaviour of the escaped references of the input vector).)]
+/
alias funRnToR1_t = real function(in real[]);

/+ FOR LATER USE.
/++
Type alias for general modifier functions which operate on the points of the
given simplex. (e.g. convert them to the nearest points in the underlying integer lattice)
+/
alias modFun_t = void delegate(ref real[][]); // Is 'ref' needed? Not sure.
+/

// ########## ENUMS ##########

/++
Type definition to mark column-major and row-major orderings of matrices. (So
far it is only used by `eigenValues()` which is a LAPACK wrapper function.
LAPACK mainly works with column-major ordering while $(B doll) works with
row-major matrices.)
+/
enum MatrixLayout
{
    /// Indicator of a column-major matrix
    ColMajor,
    /// Indicator of a row-major matrix
    RowMajor
}

// ########## STRUCTS ##########

/++
Struct for safe bound representation. The constructor performs the length
comparison on the two given vectors to guarantee the "dimensional correctness"
(i.e. the two vectors are in the same $(B n)-dimensional Euclidean space). In
the case when the two given vectors do not have the same length, the constructor
leaves the bounds empty. In this case we can change the values of the bound
vectors only with the copy (postblit) constructor. It is not allowed to change
the vector values individually. Though, if this property is still needed (and
also in the most cases if it is possible), we recommend to use the compile-
time checked version of `Bounds` (`BoundsCT`) with a given `dim` parameter.
+/
struct Bounds
{
    private
    {
        real[] _lower;
        real[] _upper;
    }

    /// Getter property for $(B lower) member (Modification is not allowed.)
    @property const(real[]) lower() const { return _lower; }

    /// Getter property for $(B upper) member (Modification is not allowed.)
    @property const(real[]) upper() const { return _upper; }

    /// Empty getter property
    @property bool empty() { return _lower.length == 0 && _upper.length == 0; }

    /// Get the dimension of the space we are working in.
    @property size_t dim() { return _lower.length; }

    /// Default constructor dupping the two slices
    /// to avoid accidental modification.
    this(real[] lower_, real[] upper_)
    {
        if (lower_.length == upper_.length)
        {
            _lower = lower_.dup;
            _upper = upper_.dup;
        }
    }

    /+/// For compile-time assert on vector lengths.
    this(size_t dim)(real[dim] lower_, real[dim] upper_)
    {
        lower = lower_;
        upper = upper_;
    }+/

    /// Postblit for deep copy.
    this(this)
    {
        _lower = _lower.dup;
        _upper = _upper.dup;
    }

    // OPERATORS

    /// Assignment operator overloading for easier
    /// initialization with arrays of $(B slices).
    void opAssign(real[][2] arr)
    {
        if (arr[0].length == arr[1].length)
        {
            _lower = arr[0].dup;
            _upper = arr[1].dup;
        }
    }

    /// Assignment operator overloading for easier
    /// initialization with arrays of $(B fixed length arrays).
    void opAssign(size_t fixDim)(real[fixDim][2] arr)
    {
        _lower.length =
        _upper.length = fixDim;

        _lower = arr[0].dup;
        _upper = arr[1].dup;
    }

    /// Equality test for two `Bounds` objects.
    bool opEquals(Bounds other) const
    {
        if (other.lower == _lower && other.upper == _upper)
            return true;
        else
            return false;
    }

    /// To ensure the two vectors share the same length
    /// (redundant --> always will be true but it is informative)
    invariant
    {
        assert(_lower.length == _upper.length);
    }
}

///
unittest
{
    Bounds bounds;
    assert(bounds.empty == true && bounds.dim == 0);

    // We can assign an array of to points to our 'Bounds' object
    bounds = [[-7, -1, 4.], [-0.2, 0, 9]];
    assert(bounds.lower == [-7, -1, 4.0L] && bounds.upper == [-0.2L, 0, 9]);
    assert(bounds.empty == false);
    assert(bounds.dim == 3);

    // Default constructor, copy and equality test
    auto copy1 = Bounds([1, 2], [3, 4]);
    auto copy2 = copy1;
    assert(copy1 == copy2);
}

/++
Compile-time checked version of Bounds.
+/
struct BoundsCT(size_t dim_)
{
    /// Lower bound vector with fixed dimension.
    real[dim_] lower;
    /// Upper bound vector with fixed dimension.
    real[dim_] upper;
    /// Get the dimension of the space we are working in.
    @property size_t dim() { return dim_; }
}

///
unittest
{
    import std.algorithm : any, map;
    import std.math : cmp;

    // We create an "empty" 'BoundsCT' object with 5-dimensional endpoints
    BoundsCT!5 bounds;
    
    // Indeed the dimension is 5
    assert(bounds.dim == 5);
    // The two vector are still uninitialized (filled with NaNs [== real.init])
    assert(!bounds.lower[].map!(a => cmp(a, real.init)).any);
    assert(!bounds.upper[].map!(a => cmp(a, real.init)).any);

    // Unlike the endpoints of a 'Bounds' object these can be modified
    bounds.lower = [1, 2, 3, 4, 5];
    bounds.upper[] = 1;

    // Now the new values should be present
    real[5] sum;
    sum[] = bounds.lower[] + bounds.upper[];
    
    assert(sum == [2.0L, 3, 4, 5, 6]);
}

/++
Struct to represent a $(LINK2 https://www.wikiwand.com/en/Simplex, proper
mathematical simplex) of $(B n+1) points in the $(B n)-dimensional space.
+/
struct Simplex
{
    import std.conv : to;

    private
    {
        real[][] _points;
        size_t _dim;
        real[] _pointSum;
    }

    // GETTER AND SETTER PROPERTIES

    /// Setter property for changing the represented point set. (Slice-version
    /// mainly for run-time calls.)
    @property void points(real[][] points_)
    {
        setPoints(points_);
    }

    /// Setter property for changing the represented point set. (Array-version
    /// mainly for compile-time calls.)
    @property void points(size_t m, size_t n)(real[n][m] points_)
    {
        setPoints(to!(real[][])(points_));
    }

    /// Getter property to query the point set defining the simplex. Be careful,
    /// it can also be used to modify the points with "slicing + assignment".
    /// This feature is intentional.
    @property real[][] points()
    {
        return _points;
    }

    /// Getter property to get the dimension of the simplex points. It is not
    /// allowed to modify this struct member explicitly only through the
    /// `randomInit()` method wich regenerates all of the simple points
    /// according to the given dimension.
    @property size_t dim()
    {
        return _dim;
    }

    /// Property for getting the sum of all points defining the simplex. It is
    /// called only in few cases from outside the struct.
    @property real[] pointSum()
    {
        calculatePointSum;
        return _pointSum.dup;
    }

    // CONSTRUCTORS

    /// If only the dimension of the simplex is given, a randomly initialized
    /// $(B n)-dimensional `Simplex` object will be returned.
    this(size_t dim_)
    {
        randomInit(dim_);
    }

    /// Postblit for deep copy.
    this(this)
    {
        _points = _points.dup;
    }

    /// Constructor for dynamic arrays destined for
    /// runtime calls (user input, data file etc.).
    this(real[][] points_)
    {
        setPoints(points_);
    }

    /// Constructor for fixed-size array inputs.
    this(size_t m, size_t n)(real[n][m] points_)
    if (m == n + 1)
    {
        setPoints(to!(real[][])(points_));
    }

    // PUBLIC AND PRIVATE MEMBER FUNCTIONS

    /// Generates $(B dim_+1) randomly initialized $(B n)-dimensional points.
    /// (Now with fixed boundaries for each coordinate: [-10, 10])
    void randomInit(size_t dim_)
    // TODO: This algorithm may be replaced with Latin Hypercube Sampling
    {
        _dim = dim_;
        _points.length = _dim + 1;
        foreach (ref pt; _points)
            pt.length = _dim;

        import std.random : uniform;

        foreach (ref coord; _points[0])
            coord = uniform(-5.0L, 5.0L);
        immutable randStep = uniform(-5.0L, 5.0L);

        foreach (direction, ref pt; _points[1 .. $])
        {
            pt = _points[0].dup;
            pt[direction] += randStep;
        }
    }

    /// Calculates the centroid (geometric center) of the point defining the
    /// simplex excluding the point with the index given as argument.
    real[] centroid(size_t excludedPointInd)
    {
        calculatePointSum;

        real[] cen;
        cen.length = _dim;

        cen[] = _pointSum[] - _points[excludedPointInd][];
        cen[] /= _dim;

        return cen;
    }

    private void setPoints(real[][] points_)
    {
        _points.length = points_.length;
        _dim = points_[0].length;

        foreach (ind, ref pt; _points)
        {
            pt.length = _dim;
            pt = points_[ind].dup;
        }
    }

    private void calculatePointSum()
    {
        _pointSum.length = _dim;
        _pointSum[] = 0.0;

        foreach (pt; _points)
            _pointSum[] += pt[];
    }

    // OPERATORS

    /// Equality test for two `Simplex` objects.
    bool opEquals(Simplex other) const
    {
        if (other.points == _points)
            return true;
        else
            return false;
    }

    /// To ensure that our `Simplex` is correct in mathematical sense.
    invariant
    {
        import std.algorithm : all;

        assert(_points.all!(a => a.length + 1 == _points.length));
    }
}

///
unittest
{
    real[][] initPoints = [[1, 2], [-2, 3], [2, 4]];
    real[2][3] fixedInit = [[1, 2], [-2, 3], [2, 4]];

    auto simplex = Simplex(initPoints);
    auto simplexFixed = Simplex(fixedInit);

    // Test the `opEquals` operator
    assert(simplex == simplexFixed);
    assert(Simplex([[]]) == Simplex([[]]));

    assert(simplex.pointSum == [1.0L, 9.0L]);

    // Centroid calculations for every index
    assert(simplex.centroid(0) == [0.0L, 3.5L]);
    assert(simplex.centroid(1) == [1.5L, 3.0L]);
    assert(simplex.centroid(2) == [-0.5L, 2.5L]);

    // Modifying points with "slicing + assignment" is allowed
    simplex.points[1][1] = -4.0L;

    // It is important to note that the original points remained unchanged
    // because of the explicit copy in the constructor
    assert(simplex.points != initPoints &&
        simplex.points == cast(real[][]) [[1, 2], [-2, -4], [2, 4]]
    );

    assert(simplex.pointSum == [1.0L, 2.0L]);
}

// ########## FUNCTIONS ##########

/// Euclidean norm for vectors given as `real` slices.
real euclideanNorm(real[] x)
{
    import std.algorithm : map, sum;
    import std.math : sqrt;
    return x.map!"a * a".sum.sqrt;
}

/++
This function generates a random Latin Hypercube Design. The dimension of
the actual space is the dimension of the given `bounds` and `numPoints`
random point will be generated. The function creates a bunch of LHDs with
the given parameters and it choose the one with the highest minimum distance
between its points. (Maximin distance)
+/
real[][] latinHypercubeDesignMaxiMin
    (Bounds bounds, size_t numPoints, size_t populationCardinality = 53)
// TODO: 'MaxiMin', 'MinCorr' versions based on a template parameter
{
    import std.algorithm : map, maxIndex, minElement;
    import std.array : array;
    import std.math : abs;
    import std.random : randomShuffle, uniform01;
    import std.range : iota;

    immutable dim = bounds.dim;

    real[][][] LatinHypercubes;
    real[] minDistances;
    LatinHypercubes.length =
    minDistances.length    = populationCardinality;

    // Store only the upper triangular matrix elements (above the diagonal)
    real[] distances;
    distances.length = numPoints * (numPoints - 1) / 2;
    size_t distInd;

    auto perm = iota(.0L, numPoints).array; // real[] (because of .0L)

    foreach (LHInd, ref LH; LatinHypercubes)
    {
        distances[] = .0L;

        // We work with the transposed LH (dim x numPoints)
        // for efficiency reasons
        LH.length = dim;
        foreach (ref dimVec; LH)
        {
            perm.randomShuffle;
            dimVec = perm.map!(a => a + uniform01).array;
            dimVec[] /= numPoints;

            // Calculate the distance respect to the current dimension
            distInd = 0;
            foreach (pointInd1, coordVal1; dimVec[0 .. $ - 1])
                foreach (coordVal2; dimVec[pointInd1 + 1 .. $])
                    // Absolute value instead of Euclidean norm
                    distances[distInd++] += abs(coordVal1 - coordVal2);
        }
        minDistances[LHInd] = distances.minElement;
    }

    auto maxiMinInd = minDistances.maxIndex;

    // Store the transposed, scaled (scalingFactor[]) and translated
    // (bounds.lower[]) Latin Hypercube in 'bestLH' and return it
    real[][] bestLH;
    bestLH.length = numPoints;

    real[] scalingFactor;
    scalingFactor.length = dim;
    scalingFactor[] = bounds.upper[] - bounds.lower[];

    foreach (pointInd, ref point; bestLH)
    {
        point.length = dim;
        foreach (dimInd, dimVec; LatinHypercubes[maxiMinInd])
            point[dimInd] = dimVec[pointInd];

        point[] *= scalingFactor[];
        point[] += bounds.lower[];
    }

    return bestLH;
}

/++
The generated LHD must contain points which are distributed among the small
hypercube parts of the bigger hypercube defined by the boundary values. Each
small hypercube should contain at most one point, and in each dimension
exactly one point should fall into each equally-sized interval between the
lowest and highest possible values of the given dimension.
+/
unittest
{
    import std.algorithm : sort;
    import std.array : array;

    immutable dim = 3;
    immutable numPoints = 4;
    auto bounds = Bounds([3, 1, -2], [7, 5, 2]); // 'dim'-dimensional bounds

    // Generate the LHD
    auto LHD = latinHypercubeDesignMaxiMin(bounds, numPoints);
    
    real[dim] distanceVec;
    distanceVec[] = bounds.upper[] - bounds.lower[];

    foreach (ind; 0 .. dim)
    {
        real[] coordinates;

        foreach (pt; LHD)
            coordinates ~= pt[ind];

        coordinates = coordinates.sort.array;

        foreach (crdInd, crd; coordinates)
        {
            assert(
                bounds.lower[ind] + crdInd * distanceVec[ind] / numPoints
                <= crd
            );
            assert(
                crd <
                bounds.lower[ind] + (crdInd + 1) * distanceVec[ind] / numPoints
            );
        }
    }
}

version (WithLapack)
{
    /++
    This is a wrapper of the LAPACK's gesvd function but later it should be replaced
    with a custom SVD implementation. Some useful resources can be found here:
    $(UL
        $(LI doll/literature/Algorithms_for_SVD.pdf)
        $(LI https://research.fb.com/fast-randomized-svd/)
        $(LI Numerical Recipes) 
    )
    +/
    double[] eigenValues(double[][] matrix, MatrixLayout layout = MatrixLayout.RowMajor)
    // TODO: Write the custom implementation! (It would be nice to implement it with
    // 'real' types istead of 'double'.)
    {
        import lapack : gesvd_, lapackint;
        import std.algorithm : joiner;
        import std.array : array;
    
        // Dimensions of the column major matrix (it may differ from the input)
        lapackint m; // Number of rows
        lapackint n; // Mober of columns
    
        // Flatten the given matrix (and transpose if needed because LAPACK works
        // with column major matrices)
        double[] flattened;
        if (MatrixLayout.RowMajor == layout)
        {
            m = cast(lapackint) matrix.length;
            n = cast(lapackint) matrix[0].length;
            flattened.length = cast(size_t) m * n;
            foreach (ind, ref elem; flattened)
                elem = matrix[ind % m][cast(size_t) ind / m];
        }
        else
        {
            m = cast(lapackint) matrix[0].length;
            n = cast(lapackint) matrix.length;
            flattened = matrix.joiner.array;
        }
        
        double* matrixPtr = flattened.ptr;
    
        // Initializing the parameters for the calculations
        lapackint lda = m;
        lapackint lwork = -1; // For workspace query
        lapackint info = void;
        double workSize;
    
        // (Init. cont.) We do not need the U and V^T matrices from SVD
        char jobChar = 'N';
        lapackint ldu = 1;
        lapackint ldvt = 1;
    
        // Array to store the eigenvalues
        double[] sigma;
        sigma.length = cast(size_t) m < n ? m : n;
        double* sigmaPtr = sigma.ptr;
    
        // Workspace query (with 'lwork == -1')
        gesvd_(jobChar, jobChar, m, n, null, lda, null, null, ldu, null, ldvt,
            &workSize, lwork, info);
    
        // Setting up the workspace based on the previous query (result in 'workSize')
        lwork = cast(lapackint) workSize;
        double[] workArr;
        workArr.length = lwork;
        double* work = workArr.ptr;
    
        // Calculating the SVD (without U and V^T matrices)
        gesvd_(jobChar, jobChar, m, n, matrixPtr, lda, sigmaPtr, null, ldu, null,
            ldvt, work, lwork, info);
    
        // Check the convergence (TODO: Handle it another way! 'sigma = [0.0]' or
        // something like this to avoid the program interruption with an exception)
        import std.exception : enforce;
        enforce(info == 0, "SVD did not converge.");
    
        return sigma;
    }
    
    // For later use
    //real[] eigenValues(size_t m, size_t n)(real[n][m] matrix){}
    
    /// We can give the input "matrix" in row major and also in column major form.
    unittest
    {
        double[][] matRM = [[2,3,4], [6,7,1], [4,2,0], [0,1,1], [1,0,8]];
        double[][] matCM = [[2,6,4,0,1], [3,7,2,1,0], [4,1,0,1,8]];
        assert(eigenValues(matRM) == eigenValues(matCM, MatrixLayout.ColMajor));
    }
}
