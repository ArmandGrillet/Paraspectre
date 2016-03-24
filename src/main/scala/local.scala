import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import java.io.File

object Local {
    // Parameters.
    val k = 7 // K'th neighbor used in local scaling.
    val minClusters = 2 // Minimal number of clusters in the dataset.
    val maxClusters = 6 // Maximal number of clusters in the dataset.

    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/0.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV.
        var matrix = breeze.linalg.csvread(matrixFile)

        // Centralizing and scale the data.
        val meancols = mean(matrix(::, *))
        // Waiting for fix scalanlp/breeze#450
        matrix = (matrix.t(::, *) - meancols.t).t
        matrix /= max(abs(matrix))

        // Build locally scaled affinity matrix.
        val distances = euclideanDistance(matrix) // Euclidean distance.
        val locScale = localScale(distances, k)
        val locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        // Build the normalized affinity matrix.
        val diagonalMatrix = sum(locallyScaledA(*, ::))
        val normalizedA = diag(pow(diagonalMatrix, -0.5)) * locallyScaledA * diag(pow(diagonalMatrix, -0.5))
        println(normalizedA)
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
                affinityMatrix(mainRow, secondRow) = -scala.math.pow(distanceMatrix(mainRow, secondRow), 2)
                affinityMatrix(mainRow, secondRow) = affinityMatrix(mainRow, secondRow) / (localScale(mainRow) * localScale(secondRow))
                affinityMatrix(mainRow, secondRow) = scala.math.exp(affinityMatrix(mainRow, secondRow))
                affinityMatrix(secondRow, mainRow) = affinityMatrix(mainRow, secondRow)
            }
        }

        return affinityMatrix
    }

    def largestPossibleGroupNumber(affinityMatrix: DenseMatrix[Double], minClusters: Int, maxClusters: Int): Int = {
        if (minClusters > maxClusters || minClusters < 2) {
            return 0
        }

        // Compute the Laplacian
        // TODO : Use CSCMatrix if sparse affinity matrix.
        val ones = DenseMatrix.ones[Double](affinityMatrix.rows, affinityMatrix.cols)
        val sumAffinityMatrix = sum(affinityMatrix(::, *))
        val diagonal = sqrt(ones :/ sumAffinityMatrix)
        val laplacian = diagonal * affinityMatrix * diagonal

        // Compute eigenvectors
        val svd.SVD(_, _, rightSingularVectors) = svd(laplacian)
        var eigenvectors = rightSingularVectors(::, 0 until maxClusters)

        // Compute eigenvalues
        var eigenvalues = diag(eigenvectors)
        eigenvalues = eigenvalues(0 until maxClusters)

        // Rotate eigenvectors
        eigenvectors = eigenvectors(::, 0 until minClusters)
        // TODO: compute the gradient of the eigenvectors alignment quality

        return 1
    }
}
