import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import breeze.plot._
import java.awt.{Color, Paint}
import scala.util.control._
import org.jfree.chart.axis.{NumberTickUnit, TickUnits}

package object serial {
    // Miscellanous functions

    def vertStack(matrix: DenseMatrix[Double], iterations: Int): DenseMatrix[Double] = {
        var stack = matrix
        var i = 0
        while (i < iterations - 1) {
            stack = DenseMatrix.vertcat(stack, matrix)
            i += 1
        }
        return stack
    }

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
        if (k > distanceMatrix.cols - 1) {
            return max(distanceMatrix(*, ::)) // Maximum distance.
        } else {
            var localScale = DenseVector.zeros[Double](distanceMatrix.cols)
            var sortedVector = IndexedSeq(0.0)

            (0 until distanceMatrix.cols).map{col =>
                sortedVector = distanceMatrix(::, col).toArray.sorted
                localScale(col) = sortedVector(k) // Kth nearest distance., the 0th neighbor is always 0 and sortedVector(1) is the first neighbor
            }

            return localScale
        }
    }

    def locallyScaledAffinityMatrix(distanceMatrix: DenseMatrix[Double], localScale: DenseVector[Double]): DenseMatrix[Double] = {
        var affinityMatrix = DenseMatrix.zeros[Double](distanceMatrix.rows, distanceMatrix.cols) // Distance matrix, size rows x cols.

        (0 until distanceMatrix.rows).map{ row =>
            (row + 1 until distanceMatrix.cols).map{ col =>
                affinityMatrix(row, col) = -scala.math.pow(distanceMatrix(row, col), 2) // -d(si, sj)²
                affinityMatrix(row, col) /= (localScale(row) * localScale(col)) // -d(si, sj)² / lambi * lambj
                affinityMatrix(row, col) = scala.math.exp(affinityMatrix(row, col)) // exp(-d(si, sj)² / lambi * lambj)
                affinityMatrix(col, row) = affinityMatrix(row, col)
            }
        }

        return affinityMatrix
    }

    // Step 5 of the self-tuning spectral clustering algorithm.

    var dims = 0
    var data = 0
    var angles = 0
    var ik, jk = DenseVector.zeros[Int](0)
    var ev = DenseMatrix.zeros[Double](0, 0)

    def paraspectre(eigenvectors: DenseMatrix[Double]): (Double, DenseVector[Int], DenseMatrix[Double]) = {
        dims = eigenvectors.cols
        data = eigenvectors.rows
        angles = (dims * (dims - 1) / 2).toInt
        ik = DenseVector.zeros[Int](angles)
        jk = DenseVector.zeros[Int](angles)
        ev = eigenvectors

        var i, j, k = 0
        for (i <- 0 until dims) {
            for (j <- (i + 1) until dims) {
                ik(k) = i
                jk(k) = j
                k += 1
            }
        }

        val maxIterations = 200

        var nablaJ, quality = 0.0
        var newQuality, old1Quality, old2Quality = 0.0 // Variables to compute the descend through true derivative.
        var qualityUp, qualityDown = 0.0 // Variables to descend through numerical derivative.
        var iter, d = 0
        var evRot = DenseMatrix.zeros[Double](0, 0)

        var theta, thetaNew = DenseVector.zeros[Double](angles)
        val loop = new Breaks

        quality = evaluateQuality(ev)
        old1Quality = quality
        old2Quality = quality

        loop.breakable{
            for (iter <- 1 to maxIterations) {
                for (d <- 0 until angles) {
                    val alpha = 0.1
                    // move up
                    thetaNew(d) = theta(d) + alpha
                    evRot = rotateGivens(thetaNew)
                    qualityUp = evaluateQuality(evRot)

                    // move down
                    thetaNew(d) = theta(d) - alpha
                    evRot = rotateGivens(thetaNew)
                    qualityDown = evaluateQuality(evRot)

                    // update only if at least one of them is better
                    if( qualityUp > quality || qualityDown > quality){
                        if( qualityUp > qualityDown ){
                            theta(d) = theta(d) + alpha
                            thetaNew(d) = theta(d)
                            quality = qualityUp
                        } else {
                            theta(d) = theta(d) - alpha
                            thetaNew(d) = theta(d)
                            quality = qualityDown
                        }
                    }
                    // val alpha = 1.0
                    // nablaJ = evaluateQualityGradient(theta, d)
                    // thetaNew(d) = theta(d) - alpha * nablaJ
                    // val evRot = rotateGivens(thetaNew)
                    // newQuality = evaluateQuality(evRot)
                    //
                    // if (newQuality > quality) {
                    //     theta(d) = thetaNew(d)
                    //     quality = newQuality
                    // } else {
                    //     thetaNew(d) = theta(d)
                    // }
                }

                if (iter > 2 && ((quality - old2Quality) < 0.001)) {
                    loop.break
                }
                old2Quality = old1Quality
                old1Quality = quality
            }
        }

        val finalEvRot = rotateGivens(thetaNew)

        if (quality equals Double.NaN) {
            return (0, null, finalEvRot)
        } else {
            val clusts = clusters(finalEvRot)
            return (quality, clusts, finalEvRot)
        }
    }

    def clusters(rotatedEigenvectors: DenseMatrix[Double]): DenseVector[Int] = {
        val absEigenvectors = abs(rotatedEigenvectors)
        return argmax(absEigenvectors(*, ::))
    }

    def evaluateQuality(matrix: DenseMatrix[Double]): Double = {
        // Take the square of all entries and find the max of each row
        var squareMatrix = pow(matrix, 2)
        val maxValues = max(squareMatrix(*, ::)) // Max of each row

        // Compute cost
        var i, j = 0
        var cost = 0.0
        for (i <- 0 until data) {
            for (j <- 0 until dims) {
                cost += squareMatrix(i, j) / maxValues(i)
            }
        }

        return 1.0 - (cost / data - 1.0) / dims
    }

    def evaluateQualityGradient(theta: DenseVector[Double], angle: Int): Double = {
        // Build V, U, A
        var vForAngle = DenseMatrix.zeros[Double](dims, dims)
        vForAngle(ik(angle),ik(angle)) = -sin(theta(angle))
        vForAngle(ik(angle),jk(angle)) = cos(theta(angle))
        vForAngle(jk(angle),ik(angle)) = -cos(theta(angle))
        vForAngle(jk(angle),jk(angle)) = -sin(theta(angle))
        val u1 = uAB(theta, 1, angle - 1)
        val u2 = uAB(theta, angle + 1, angles -1)

        val a = ev * u1 * vForAngle * u2

        val y = rotateGivens(theta)

        val maxValues = max(y(*, ::)) // Max of each row
        val maxIndexCol = argmax(y(*, ::))

        // Compute gradient
        var nablaJ, tmp1, tmp2 = 0.0
        var i, j = 0
        for (i <- 0 until data) { // Loop over all rows
            for (j <- 0 until dims) { // Loop over all columns
                tmp1 = a(i, j) * y(i, j) / (maxValues(i) * maxValues(i))
                tmp2 = a(i, maxIndexCol(i)) * pow(y(i, j), 2) / pow(maxValues(i), 3)
                nablaJ += tmp1 - tmp2
            }
        }
        nablaJ = 2 * nablaJ / data / dims

        return nablaJ
    }

    def rotateGivens(theta: DenseVector[Double]): DenseMatrix[Double] = {
        val g = uAB(theta, 0, angles - 1)
        return ev * g
    }

    def uAB(theta: DenseVector[Double], a: Int, b: Int): DenseMatrix[Double] = {
        var i, k = 0
        var uab = DenseMatrix.eye[Double](dims)

        if (b < a) {
            return uab
        }

        var tt, uIk = 0.0
        for (k <- a to b) {
            tt = theta(k)
            for (i <- 0 until dims) {
                uIk = uab(i, ik(k)) * cos(tt) - uab(i, jk(k)) * sin(tt)
                uab(i, jk(k)) = uab(i, ik(k)) * sin(tt) + uab(i, jk(k)) * cos(tt)
                uab(i, ik(k)) = uIk
            }
        }

        return uab
    }
}
