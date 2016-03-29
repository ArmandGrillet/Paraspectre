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
        var localScale = DenseVector.zeros[Double](distanceMatrix.cols)
        var sortedVector = IndexedSeq(0.0)

        if (k >= distanceMatrix.cols - 1) {
            (0 until distanceMatrix.cols).map{col =>
                localScale(col) = distanceMatrix(::, col).max // Maximum distance.
            }
        } else {
            (0 until distanceMatrix.cols).map{col =>
                sortedVector = distanceMatrix(::, col).toArray.sorted
                localScale(col) = sortedVector(k) // Kth nearest distance.
            }
        }

        return localScale
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

    def cbest(eigenvectors: DenseMatrix[Double]) : {val clusterAssignment : DenseMatrix[Double]; val quality : Double; val rotatedEigenvectors : DenseMatrix[Double]} = {
        new {
            val clusterAssignment = DenseMatrix.zeros[Double](2, 2)
            val quality = 1.0;
            val rotatedEigenvectors = DenseMatrix.zeros[Double](2, 2)
        }
    }
}
