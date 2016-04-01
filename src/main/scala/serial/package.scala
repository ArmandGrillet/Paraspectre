import breeze.linalg._
import breeze.numerics._
import breeze.stats._

package object serial {
    def euclideanDistance(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var distanceMatrix = DenseMatrix.zeros[Double](matrix.rows, matrix.rows) // Distance matrix, size rows x rows.
        var distanceVector = DenseVector(0.0).t // The distance vector containing the distance between two vectors.

        (0 until matrix.rows).map{ mainRow =>
            (mainRow + 1 until matrix.rows).map{ secondRow =>
                distanceVector = matrix(mainRow, ::) - matrix(secondRow,::) // Xi - Xj | Yi - Yj
                distanceVector *= distanceVector // (Xi - Xj)² | (Yi - Yj)²
                distanceMatrix(mainRow, secondRow) = sqrt(sum(distanceVector)) // √(Xi - Xj)² + (Yi - Yj)² + ...
                distanceMatrix(secondRow, mainRow) = distanceMatrix(mainRow, secondRow)
            }
        }

        return distanceMatrix
    }

    def localScale(distanceMatrix: DenseMatrix[Double], k: Int): DenseVector[Double] = {
        if (k >= distanceMatrix.cols - 1) {
            return max(distanceMatrix(*, ::)) // Maximum distance.
        } else {
            var localScale = DenseVector.zeros[Double](distanceMatrix.cols)
            var sortedVector = IndexedSeq(0.0)

            (0 until distanceMatrix.cols).map{col =>
                sortedVector = distanceMatrix(::, col).toArray.sorted
                localScale(col) = sortedVector(k) // Kth nearest distance.
            }

            return localScale
        }
    }

    def locallyScaledAffinityMatrix(distanceMatrix: DenseMatrix[Double], localScale: DenseVector[Double]): DenseMatrix[Double] = {
        var affinityMatrix = DenseMatrix.zeros[Double](distanceMatrix.rows, distanceMatrix.rows) // Distance matrix, size rows x rows.

        (0 until distanceMatrix.rows).map{ mainRow =>
            (mainRow + 1 until distanceMatrix.rows).map{ secondRow =>
                affinityMatrix(mainRow, secondRow) = -scala.math.pow(distanceMatrix(mainRow, secondRow), 2) // -d(si, sj)²
                affinityMatrix(mainRow, secondRow) = affinityMatrix(mainRow, secondRow) / (localScale(mainRow) * localScale(secondRow)) // -d(si, sj)² / lambi * lambj
                affinityMatrix(mainRow, secondRow) = scala.math.exp(affinityMatrix(mainRow, secondRow)) // exp(-d(si, sj)² / lambi * lambj)
                affinityMatrix(secondRow, mainRow) = affinityMatrix(mainRow, secondRow)
            }
        }

        return affinityMatrix
    }

    def logicalNot(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var result = DenseMatrix.zeros[Double](matrix.rows, matrix.cols)
        var row = 0
        while (row < matrix.rows) {
            var col = 0
            while (col < matrix.cols) {
                if (matrix(row, col) == 0.0) {
                    result(row, col) = 1.0
                }
                col += 1
            }
            row += 1
        }

        return result
    }

    def rotateEigenvectors(eigenvectors: DenseMatrix[Double]): (DenseMatrix[Double], Double, DenseMatrix[Double]) = {
        // Only works in two dimensions.
        val ndata = eigenvectors.rows
        val dim = eigenvectors.cols

        // Get the number of angles
        val angles = (dim*(dim-1)/2).toInt
        val theta = DenseVector.zeros[Double](angles)

        // We know that the method is 1

        // Build index mapping
        var ik = DenseVector.zeros[Int](angles)
        var jk = DenseVector.zeros[Int](angles)
        var k = 0
        var i = 0
        while (i < dim - 1) {
            var j = i + 1
            while (j < dim) {
                ik(k) = i
                jk(k) = j
                k += 1
                j += 1
            }
            i += 1
        }

        // Definitions
        val maxIterations = 200
        var iteration, d = 0
        var q, qOld1, qOld2 = alignmentQuality(eigenvectors)

        while (iteration < maxIterations) {
            iteration += 1
            while (d < angles) {
                val dQ = qualityGradient(eigenvectors, theta, ik, jk, angles, d)
                d += 1
            }
        }


        return (DenseMatrix.zeros[Double](2, 2), 1.0, DenseMatrix.zeros[Double](2, 2))
    }

    def alignmentQuality(eigenvectors: DenseMatrix[Double]): Double = {
        // Take the square of all entries and find the max of each row
        val squaredMatrix = eigenvectors :* eigenvectors // :* = Hadamard product
        val maxEachRow = max(squaredMatrix(*, ::))

        // Compute cost
        var cost = 0.0
        var row = 0
        while (row < eigenvectors.rows) {
            var col = 0
            while (col < eigenvectors.cols) {
                cost += squaredMatrix(row, col) / maxEachRow(row)
                col += 1
            }
            row += 1
        }

        cost = 1.0 - (cost / eigenvectors.rows - 1.0) / eigenvectors.cols

        return cost
    }

    def qualityGradient(eigenvectors: DenseMatrix[Double], theta: DenseVector[Double], ik: DenseVector[Int], jk: DenseVector[Int], angles: Int, index: Int): Double = {
        // In C++ it is a 1D array but arr(rows * row + col) = mat(col, row)
        var gradients = DenseMatrix.zeros[Double](eigenvectors.cols, eigenvectors.cols)
        gradients(ik(index), ik(index)) = -sin(theta(index))
        gradients(jk(index), ik(index)) = cos(theta(index))
        gradients(ik(index), jk(index)) = -cos(theta(index))
        gradients(jk(index), jk(index)) = -sin(theta(index))

        val u1 = uAB(theta, 0, index - 1, ik, jk, eigenvectors.cols)
        val u2 = uAB(theta, index + 1, angles - 1, ik, jk, eigenvectors.cols)

        val a = eigenvectors * gradients * u1 * u2
        return 1.0
    }

    def uAB(theta: DenseVector[Double], a: Int, b: Int, ik: DenseVector[Int], jk: DenseVector[Int], dim: Int): DenseMatrix[Double] = {
        // Set uab to be an identity matrix
        var uab = DenseMatrix.eye[Double](dim)

        if (b < a) {
            return uab
        }

        var i = 0
        var k = a
        var tt, uIndex = 0.0
        while (k <= b) {
            tt = theta(k)
            while (i < dim) {
                uIndex = uab(i, ik(k)) * cos(tt) - uab(i, jk(k)) * sin(tt)
                uab(i, jk(k)) = uab(i, ik(k)) * sin(tt) * uab(i, jk(k)) * cos(tt)
                uab(i, ik(k)) = uIndex
                i += 1
            }
            k += 1
        }

        return uab
    }
}
