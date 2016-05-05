package serial

import breeze.linalg._
import breeze.numerics._
import breeze.plot._
import breeze.stats._
import java.awt.{Color, Paint}
import java.io.File
import scala.collection.mutable.ListBuffer
import scala.io.Source

class Algorithm(argDataset: DenseMatrix[Double], argMinClusters: Int, argMaxClusters: Int, argDebug: Boolean) {
    // Parameters.
    val k = 7 // Kth neighbor used in local scaling.
    val dataset = argDataset
    val minClusters = argMinClusters // Minimal number of clusters in the dataset.
    val maxClusters = argMaxClusters // Maximal number of clusters in the dataset.
    val printer = new Printer(argDebug)

    def cluster(): DenseVector[Int] = {

        // Centralize and scale the data.
        // val meanCols = mean(originalMatrix(::, *)).t.toDenseMatrix
        // var matrix = (originalMatrix - vertStack(meanCols, originalMatrix.rows))
        // matrix /= max(abs(matrix))
        //val matrix = originalMatrix

        printer.printDataset(dataset)

        // Compute local scale (step 1).
        val distances = euclideanDistance(dataset)
        val locScale = localScale(distances, k)

        // Build locally scaled affinity matrix (step 2).
        var locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        // Build the normalized affinity matrix (step 3)
        val diagonalMatrix = diag(pow(sum(locallyScaledA(*, ::)), -0.5)) // Sum of each row, then power -0.5, then matrix.
        val normalizedA = diagonalMatrix * locallyScaledA * diagonalMatrix

        // Compute the largest eigenvectors
        val eigenstuff = eig(normalizedA)
        var eigenvalues = eigenstuff.eigenvalues // DenseVector
        val unsortedEigenvectors = eigenstuff.eigenvectors // DenseMatrix
        var eigenvectors = DenseMatrix.zeros[Double](unsortedEigenvectors.rows, maxClusters)
        var biggestEigenvalues = DenseVector.zeros[Double](maxClusters)

        var i = 0
        val minEigenvalue = min(eigenvalues)
        for (i <- 0 until maxClusters) {
            val indexBiggestEigenvalue = argmax(eigenvalues)
            biggestEigenvalues(i) = eigenvalues(indexBiggestEigenvalue)
            eigenvalues(indexBiggestEigenvalue) = minEigenvalue
            for (row <- 0 until unsortedEigenvectors.rows) {
                eigenvectors(row, i) = unsortedEigenvectors(row, indexBiggestEigenvalue)
            }
        }

        printer.printEigenvalues(biggestEigenvalues)

        // In cluster_rotate.m originally
        var qualities = new ListBuffer[Double]()
        var currentEigenvectors = eigenvectors(::, 0 until minClusters)
        var (quality, clusters, rotatedEigenvectors) = paraspectre(currentEigenvectors)
        qualities += quality

        println(minClusters + " clusters:\t" + quality)

        var group = 0
        for (group <- minClusters until maxClusters) {
            val eigenvectorToAdd = eigenvectors(::, group).toDenseMatrix.t
            currentEigenvectors = DenseMatrix.horzcat(rotatedEigenvectors, eigenvectorToAdd)
            val (tempQuality, tempClusters, tempRotatedEigenvectors) = paraspectre(currentEigenvectors)
            qualities += tempQuality
            rotatedEigenvectors = tempRotatedEigenvectors
            println((group + 1) + " clusters:\t" + tempQuality)

            if (tempQuality >= quality - 0.001) {
                quality = tempQuality
                clusters = tempClusters
            }
        }

        printer.printQualities(qualities, minClusters)
        printer.printClusters(dataset, clusters)
        return clusters
    }
}
