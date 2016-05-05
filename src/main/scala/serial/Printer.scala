package serial

import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import breeze.plot._
import java.awt.{Color, Paint}
import scala.collection.mutable.ListBuffer
import scala.util.control._
import org.jfree.chart.axis.{NumberTickUnit, TickUnits}

class Printer(argDebug: Boolean) {
    val debug = argDebug
    val f = Figure()

    def printDataset(dataset: DenseMatrix[Double]) {
        if (debug) {
            val p = f.subplot(0)
            p.title = "Dataset"
            p += scatter(dataset(::, 0), dataset(::, 1), {(_:Int) => 0.01}, {(pos:Int) => Color.BLACK}) // Display the observations.
        }
    }

    def printLocalScaling(dataset: DenseMatrix[Double], locallyScaledAffinityMatrix: DenseMatrix[Double]) {
        if (debug) {
            val p = f.subplot(0)
            p += scatter(dataset(::, 0), dataset(::, 1), {(_:Int) => 0.01}, {(pos:Int) => Color.RED}) // Display the observations.

            var sortedVector = IndexedSeq(0.0)
            (0 until dataset.rows).map{ firstRow =>
                sortedVector = locallyScaledAffinityMatrix(::, firstRow).toArray.sorted // I should do distanceMatrix(row, ::) but as it is a distance matrix thus it is equivalent.
                (firstRow + 1 until dataset.rows).map{ secondRow =>
                    var minRow, maxRow = 0
                    if (dataset(firstRow, 0) <= dataset(secondRow, 0)) {
                        minRow = firstRow
                        maxRow = secondRow
                    } else {
                        minRow = secondRow
                        maxRow = firstRow
                    }
                    val x = linspace(dataset(minRow, 0), dataset(maxRow, 0))
                    val factor = (dataset(maxRow, 1) - dataset(minRow, 1)) / (dataset(maxRow, 0) - dataset(minRow, 0))
                    val y = DenseVector.tabulate(x.length){i => dataset(minRow, 1) + factor * (x(i) - dataset(minRow, 0))}
                    p += plot(x, y)
                }
            }
        }
    }

    def printEigenvalues(eigenvalues: DenseVector[Double]) {
        if (debug) {
            val p = f.subplot(2,2,2)
            p.title = "First " + eigenvalues.length + " eigenvalues of L"
            p.yaxis.setTickUnit(new NumberTickUnit(0.005))
            val xVector = linspace(0, eigenvalues.length - 1, eigenvalues.length)

            p += scatter(xVector, eigenvalues, {(_:Int) => 0.05}, {(_:Int) => Color.RED}) // Display the observations.

            (0 until (eigenvalues.length - 1)).map{ eigenvalue =>
                val x = linspace(eigenvalue, eigenvalue + 1)
                val factor = eigenvalues(eigenvalue + 1) - eigenvalues(eigenvalue)
                val y = DenseVector.tabulate(x.length){i => eigenvalues(eigenvalue) + factor * (x(i) - x(0))}
                p += plot(x, y, '-', "BLUE")
            }
        }
    }

    def printClusters(matrix: DenseMatrix[Double], clusters: DenseVector[Int]) {
        if (debug) {
            val colors = List(Color.RED, Color.GREEN, Color.BLUE, Color.BLACK, Color.MAGENTA, Color.CYAN, Color.YELLOW)

            val p = f.subplot(2,2,1)
            p.title = "Clusters"

            p += scatter(matrix(::, 0), matrix(::, 1), {(_:Int) => 0.01}, {(pos:Int) => colors(clusters(pos))}) // Display the observations.
        }
    }

    def printQualities(qualities: ListBuffer[Double], minClusters: Int) {
        if (debug) {
            val p = f.subplot(2,2,3)
            p.title = "Qualities"
            p.yaxis.setTickUnit(new NumberTickUnit(0.01))

            val qualitiesVector = new DenseVector(qualities.toArray)
            val xVector = linspace(minClusters, (qualities.length - 1) + minClusters, qualities.length)

            p += scatter(xVector, qualitiesVector, {(_:Int) => 0.05}, {(_:Int) => Color.RED}) // Display the observations.

            (0 until (qualities.length - 1)).map{ quality =>
                val x = linspace(minClusters + quality, minClusters + quality + 1)
                val factor = qualities(quality + 1) - qualities(quality)
                val y = DenseVector.tabulate(x.length){i => qualities(quality) + factor * (x(i) - x(0))}
                p += plot(x, y, '-', "BLUE")
            }
        }
    }
}
