package serial

import breeze.linalg._
import breeze.numerics._
import breeze.plot._
import breeze.stats._
import java.awt.{Color, Paint}
import java.io.File
import scala.io.Source

object Algorithm {
    // Parameters.
    val k = 7 // K'th neighbor used in local scaling.
    val minClusters = 2 // Minimal number of clusters in the dataset.
    val maxClusters = 6 // Maximal number of clusters in the dataset.
    val eps = 2.2204e-16

    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/4.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV.
        val originalMatrix = breeze.linalg.csvread(matrixFile)

        // Centralizing and scale the data.
        val meanCols = mean(originalMatrix(::, *)).t.toDenseMatrix
        var matrix = (originalMatrix - vertStack(meanCols, originalMatrix.rows))
        matrix /= max(abs(matrix))

        // Build locally scaled affinity matrix (step 2).
        val distances = euclideanDistance(matrix)
        val locScale = localScale(distances, k)
        var locallyScaledA = locallyScaledAffinityMatrix(distances, locScale)

        // Build the normalized affinity matrix.
        val diagonalMatrix = sum(locallyScaledA(*, ::)) // Sum of each row (eigenvector)
        val normalizedA = diag(pow(diagonalMatrix, -0.5)) * locallyScaledA * diag(pow(diagonalMatrix, -0.5))

        // In evecs.m originally
        // Compute the Laplacian
        // TODO : Use CSCMatrix if sparse affinity matrix.
        // Sum down each column
        val sumCols = (sum(normalizedA(::, *)) :+ eps).t
        val diagonal = diag(sqrt(DenseVector.ones[Double](sumCols.length) :/ sumCols))
        val laplacian = diagonal * normalizedA * diagonal

        // Compute eigenvectors
        val svd.SVD(_, _, rightSingularVectors) = svd(laplacian)
        val eigenvectors = rightSingularVectors(::, 0 to maxClusters)

        // Compute the eigenvalues
        // var eigenvalues = diag(eigenvectors)
        // eigenvalues = eigenvalues(0 until maxClusters)
        // end of evecs.m

        // In cluster_rotate.m originally
        var currentEigenvectors = eigenvectors(::, 0 until minClusters)
        var (cost, clusters, rotatedEigenvectors) = paraspectre(currentEigenvectors)

        print(minClusters)
        print(" clusters:\t")
        println(cost)

        var group = 0
        for (group <- (minClusters + 1) to maxClusters) {
            val eigenvectorToAdd = eigenvectors(::, group).toDenseMatrix.t
            currentEigenvectors = DenseMatrix.horzcat(rotatedEigenvectors, eigenvectorToAdd)
            val (tempCost, tempClusters, tempRotatedEigenvectors) = paraspectre(currentEigenvectors)
            rotatedEigenvectors = tempRotatedEigenvectors
            print(group)
            print(" clusters:\t")
            println(tempCost)
            if (tempCost <= (cost + 0.001)) {
                cost = tempCost
                clusters = tempClusters
            }
        }
        // In evrot.cpp originally

        // val f = Figure()
        // val id2Color: Int => Paint = id => id match {
        //     case 0 => Color.YELLOW
        //     case 1 => Color.RED
        //     case 2 => Color.GREEN
        //     case 3 => Color.BLUE
        //     case 4 => Color.GRAY
        //     case _ => Color.BLACK
        //   }
        //
        // f.subplot(0) +=  scatter(originalMatrix(::, 0), originalMatrix(::, 1), {(_:Int) => 1.0}, {(_:Int) => Color.BLACK})
        // f.subplot(0).xlabel = "X-coordinate"
        // f.subplot(0).ylabel = "Y-coordinate"
    }
}
