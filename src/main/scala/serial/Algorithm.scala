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
    val eps = 2.2204e-16

    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/0.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV.
        var matrix = breeze.linalg.csvread(matrixFile)

        // Centralizing and scale the data.
        val meancols = mean(matrix(::, *))
        matrix = (matrix.t(::, *) - meancols.t).t // Waiting for fix scalanlp/breeze#450
        matrix /= max(abs(matrix))

        // Build locally scaled affinity matrix.
        val distances = euclideanDistance(matrix)
        val locScale = localScale(distances, k)
        var locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        // Zero out diagonal
        locallyScaledA = locallyScaledA :* logicalNot(DenseMatrix.eye[Double](matrix.rows)) // logicalNot = ~

        // Build the normalized affinity matrix, in the paper but not the code.
        // val diagonalMatrix = sum(locallyScaledA(*, ::))
        // val normalizedA = diag(pow(diagonalMatrix, -0.5)) * locallyScaledA * diag(pow(diagonalMatrix, -0.5))

        // In evecs.m originally
        // Compute the Laplacian
        // TODO : Use CSCMatrix if sparse affinity matrix.
        // Sum down each column
        val sumAffinityMatrix = sum(locallyScaledA(::, *)) + DenseVector.fill(locallyScaledA.rows){eps}.t
        val diagonal = diag(sqrt(DenseVector.ones[Double](locallyScaledA.rows).t :/ sumAffinityMatrix))
        val laplacian = diagonal * locallyScaledA * diagonal

        // Compute eigenvectors
        val svd.SVD(_, _, rightSingularVectors) = svd(laplacian)
        val eigenvectors = rightSingularVectors(::, 0 until maxClusters)

        // Compute the eigenvalues
        // var eigenvalues = diag(eigenvectors)
        // eigenvalues = eigenvalues(0 until maxClusters)
        // end of evecs.m

        // In cluster_rotate.m originally
        var qualities = scala.collection.mutable.MutableList[Double]()
        var currentEigenvectors = eigenvectors(::, 0 until minClusters)
        var (quality, clusters, rotatedEigenvectors) = rotateEigenvectors(currentEigenvectors)
        qualities += quality

        var group = 0
        for (group <- (minClusters + 1) to maxClusters) {
            currentEigenvectors = DenseMatrix.horzcat(rotatedEigenvectors, eigenvectors(::, 0 until group))
            val (tempQuality, tempClusters, tempRotatedEigenvectors) = rotateEigenvectors(currentEigenvectors)
            rotatedEigenvectors = tempRotatedEigenvectors
            qualities += tempQuality
        }

        val i = qualities.filter(quality => max(qualities) - quality <= 0.001)
        val bestGroupIndex = i.last
        println(bestGroupIndex)
        // In evrot.cpp originally



        // Rotate eigenvectors
        // var currentEigenvectors = eigenvectors(::, 0 until minClusters)

        // TODO: compute the gradient of the eigenvectors alignment quality
    }
}
