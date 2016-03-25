package serial

import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import java.io.File

object Algorithm {
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
}
