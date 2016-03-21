import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import scala.math.{pow, exp}
import java.io.File

object Local {
    // Parameters.
    val k = 7 // K'th neighbor used in local scaling.

    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/0.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV.
        var matrix = breeze.linalg.csvread(matrixFile)
        // println(matrix)

        // Centralizing and scale the data.
        val meancols = mean(matrix(::, *))
        // Waiting for fix scalanlp/breeze#450
        matrix = (matrix.t(::, *) - meancols.t).t
        matrix /= max(abs(matrix))

        // Build affinity matrix.
        val distances = euclideanDistance(matrix) // Euclidean distance.
        val locScale = localScale(distances, k)
        val locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        println(locallyScaledA)
    }

    def euclideanDistance(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var distanceMatrix = DenseMatrix.zeros[Double](matrix.rows, matrix.rows) // Distance matrix, size rows x rows.
        var distanceVector = DenseVector(0.0).t // The distance vector containing the distance between two vectors.

        (0 until matrix.rows).map{ mainRow =>
            (mainRow + 1 until matrix.rows).map{ secondRow =>
                distanceVector = matrix(mainRow, ::) - matrix(secondRow,::) // Xi - Xj | Yi - Yj
                distanceVector *= distanceVector // (Xi - Xj)^^2 | (Yi - Yj)^^2
                distanceMatrix(mainRow, secondRow) = sqrt(sum(distanceVector)) // √(Xi - Xj)^^2 + (Yi - Yj)^^2 + ...
                distanceMatrix(secondRow, mainRow) = distanceMatrix(mainRow, secondRow)
            }
        }

        return distanceMatrix
    }

    def localScale(distanceMatrix: DenseMatrix[Double], k: Int): DenseVector[Double] = {
        var localScale = DenseVector.zeros[Double](distanceMatrix.cols)
        var sortedVector = IndexedSeq(0.0)
        if (k > distanceMatrix.cols) {
            (0 until distanceMatrix.cols).map{col =>
                localScale(col) = distanceMatrix(::, col).max
            }
        } else {
            (0 until distanceMatrix.cols).map{col =>
                sortedVector = distanceMatrix(::, col).toArray.sorted
                localScale(col) = sortedVector(k - 1)
            }
        }

        return localScale
    }

    def locallyScaledAffinityMatrix(distanceMatrix: DenseMatrix[Double], localScale: DenseVector[Double]): DenseMatrix[Double] = {
        var affinityMatrix = DenseMatrix.zeros[Double](distanceMatrix.rows, distanceMatrix.rows) // Distance matrix, size rows x rows.

        (0 until distanceMatrix.rows).map{ mainRow =>
            (mainRow + 1 until distanceMatrix.rows).map{ secondRow =>
                affinityMatrix(mainRow, secondRow) = -pow(distanceMatrix(mainRow, secondRow), 2)
                affinityMatrix(mainRow, secondRow) = affinityMatrix(mainRow, secondRow) / (localScale(mainRow) * localScale(secondRow))
                affinityMatrix(mainRow, secondRow) = exp(affinityMatrix(mainRow, secondRow))
                affinityMatrix(secondRow, mainRow) = affinityMatrix(mainRow, secondRow)
            }
        }

        return affinityMatrix
    }
}
