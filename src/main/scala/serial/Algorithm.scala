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
        // val meanCols = mean(dataset(::, *)).t.toDenseMatrix
        // var matrix = (dataset - vertStack(meanCols, dataset.rows))
        // matrix /= max(abs(matrix))
        // val matrix = dataset

        printer.printDataset(dataset)

        // Compute local scale (step 1).
        val distances = euclideanDistance(dataset)
        val locScale = localScale(distances, k)

        // Build locally scaled affinity matrix (step 2).
        val locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        // Build the normalized affinity matrix (step 3)
        val diagonalMatrix = diag(pow(sum(locallyScaledA(*, ::)), -0.5)) // Sum of each row, then power -0.5, then matrix.
        var normalizedA = diagonalMatrix * locallyScaledA * diagonalMatrix

        // Breeze does not think normalizedA is symmetric, let's fix that.
        var row, col = 0
        for (row <- 0 until normalizedA.rows) {
            for (col <- row + 1 until normalizedA.cols) {
                normalizedA(col, row) = normalizedA(row, col)
            }
        }

        // Compute the largest eigenvectors
        val eigenstuff = eigSym(normalizedA)
        var eigenvectors = DenseMatrix.zeros[Double](eigenstuff.eigenvectors.rows, maxClusters)
        var biggestEigenvalues = DenseVector.zeros[Double](maxClusters)

        var i = 0
        for (i <- 0 until maxClusters) {
            biggestEigenvalues(i) = eigenstuff.eigenvalues(-(1 + i))
            println("eigenvalue n°" + i + " : " + biggestEigenvalues(i))
            for (row <- 0 until eigenstuff.eigenvectors.rows) {
                eigenvectors(row, i) = eigenstuff.eigenvectors(row, -(1 + i))
            }
        }
        println("")
        printer.printEigenvalues(biggestEigenvalues)

        // csvwrite(new File("eigenvectors3.txt"), eigenvectors3, separator = ',')
        // val staticEigenvectors = new File("eigenvectors3.txt")
        // eigenvectors = breeze.linalg.csvread(staticEigenvectors)

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

            if (tempQuality >= quality - 0.002) {
                quality = tempQuality
                clusters = tempClusters
            }
        }

        printer.printQualities(qualities, minClusters)
        printer.printClusters(dataset, clusters)
        return clusters
    }
}
