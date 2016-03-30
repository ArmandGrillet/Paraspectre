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
                if (matrix(row, col) == 0) {
                    result(row, col) = 1
                }
                col += 1
            }
            row += 1
        }

        return result
    }

    def rotateEigenvectors(eigenvectors: DenseMatrix[Double]): (DenseMatrix[Double], Double, DenseMatrix[Double]) = {
        // Only works in two dimensions now.
        val ndata = eigenvectors.rows
        val dim = eigenvectors.cols

        // Get the number of angles
        val angles = (dim*(dim-1)/2).toInt
        val theta = DenseMatrix.zeros[Double](1, angles)

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
        var iteration = 0
        var q, qOld1, qOld2 = alignmentQuality(eigenvectors)


        return (DenseMatrix.zeros[Double](2, 2), 1.0, DenseMatrix.zeros[Double](2, 2))
    }

    def alignmentQuality(eigenvectors: DenseMatrix[Double]): Double = {
        // Take the square of all entries and find the max of each row
        val squareMatrix = square(eigenvectors)
        val maxEachRow = max(squareMatrix(*, ::))

        // Compute cost
        var cost = 0.0
        var row = 0
        while (row < eigenvectors.rows) {
            var col = 0
            while (col < eigenvectors.cols) {
                cost += squareMatrix(row, col) / maxEachRow(row)
                col += 1
            }
            row += 1
        }

        cost = 1.0 - (cost / eigenvectors.rows - 1.0) / eigenvectors.cols

        return cost
    }

    def square(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var result = DenseMatrix.zeros[Double](matrix.rows, matrix.cols)
        var row = 0
        while (row < matrix.rows) {
            var col = 0
            while (col < matrix.cols) {
                result(row, col) = matrix(row, col) * matrix(row, col)
                col += 1
            }
            row += 1
        }

        return result
    }

}