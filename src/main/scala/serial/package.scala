import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import breeze.plot._
import java.awt.{Color, Paint}
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

    def printVector(vector: DenseVector[Double]) {
        val f = Figure()
        val p = f.subplot(0)
        p.title = "First 10 eigenvalues of L"
        p.xlim(0, vector.length - 1)
        p.ylim(0.9, 1.01)
        p.yaxis.setTickUnit(new NumberTickUnit(0.01));

        val xVector = linspace(0, vector.length - 1, vector.length)

        p += scatter(xVector, vector, {(_:Int) => 0.3}, {(_:Int) => Color.RED}) // Display the observations.
    }

    def printClusters(matrix: DenseMatrix[Double], clusters: DenseVector[Int]) {
        val colors = List(Color.RED, Color.GREEN, Color.BLUE, Color.BLACK, Color.MAGENTA, Color.CYAN, Color.YELLOW)

        val f = Figure()
        val p = f.subplot(0)
        p.title = "Clusters"

        p += scatter(matrix(::, 0), matrix(::, 1), {(_:Int) => 0.01}, {(pos:Int) => colors(clusters(pos))}) // Display the observations.
    }

    def logicalNot(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        return DenseMatrix.tabulate(matrix.rows, matrix.cols){
            case (row, col) =>
            if (matrix(row, col) == 0.0) {
                1.0
            } else {
                0.0
            }
        }
    }

    // Self-tuning spectral clustering last steps
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
        val alpha = 1.0
        var dQ, q, qNew, qOld1, qOld2, qUp, qDown = 0.0
        var iter, d = 0

        var theta, thetaNew = DenseVector.zeros[Double](angles)

        q = evaluateQuality(ev)
        qOld1 = q
        qOld2 = q
        while (iter < maxIterations) {
            iter += 1
            for (d <- 0 until angles) {
                dQ = evaluateQualityGradient(theta, d)
                thetaNew(d) = theta(d) - alpha * dQ
                val evRot = rotateGivens(thetaNew)
                qNew = evaluateQuality(evRot)

                if (qNew > q) {
                    theta(d) = thetaNew(d)
                    q = qNew
                } else {
                    thetaNew(d) = theta(d)
                }
            }

            if (iter > 2 && ((q - qOld2) < 1e-3)) {
                iter = maxIterations
            } else {
                qOld2 = qOld1
                qOld1 = q
            }
        }

        val finalEvRot = rotateGivens(thetaNew)
        val clusts = clusters(finalEvRot)
        return (q, clusts, finalEvRot)
    }

    def clusters(rotatedEigenvectors: DenseMatrix[Double]): DenseVector[Int] = {
        val squaredVectors = rotatedEigenvectors :* rotatedEigenvectors
        return argmax(squaredVectors(*, ::))
    }

    def evaluateQuality(x: DenseMatrix[Double]): Double = {
        // Take the square of all entries and find the max of each row
        var x2 = pow(x, 2)
        val maxValues = max(x2(*, ::)) // Max of each row

        // Compute cost
        var row, col = 0
        for (row <- 0 until data) {
            for (col <- 0 until dims) {
                x2(row, col) /= maxValues(row)
            }
        }

        val j = 1.0 - (sum(x2) / data - 1.0) / dims
        return j
    }

    def evaluateQualityGradient(theta: DenseVector[Double], angle: Int): Double = {
        // Build V, U, A
        val v = gradU(theta, angle)
        val u1 = uAB(theta, 0, angle - 1)
        val u2 = uAB(theta, angle + 1, angles -1)

        val a = ev * u1 * v * u2

        val y = rotateGivens(theta)

        val maxValues = max(y(*, ::)) // Max of each row
        val maxIndexCol = argmax(y(*, ::))

        // Compute gradient
        var dJ, tmp1, tmp2 = 0.0
        var i, j = 0
        for (j <- 0 until dims) { // Loop over all columns
            for (i <- 0 until data) { // Loop over all rows
                tmp1 = a(i, j) * y(i, j) / (maxValues(i) * maxValues(i))
                tmp2 = a(i, maxIndexCol(i)) * pow(y(i, j), 2) / pow(maxValues(i), 3)
                dJ += tmp1 - tmp2
            }
        }
        dJ = 2 * dJ / data / dims

        return dJ
    }

    def rotateGivens(theta: DenseVector[Double]): DenseMatrix[Double] = {
        val g = uAB(theta, 0, angles -1)
        val y = ev * g
        return y
    }

    def uAB(theta: DenseVector[Double], a: Int, b: Int): DenseMatrix[Double] = {
        var i, k = 0
        var uab = DenseMatrix.eye[Double](dims)

        if (b < a) {
            return uab
        }

        var tt, uIk = 0.0
        for (k <- 0 to b) {
            tt = theta(k)
            for (i <- 0 until dims) {
                uIk = uab(i, ik(k)) * cos(tt) - uab(i, jk(k)) * sin(tt)
                uab(i, jk(k)) = uab(i, ik(k)) * sin(tt) + uab(i, jk(k)) * cos(tt)
                uab(i, ik(k)) = uIk
            }
        }

        return uab
    }

    def gradU(theta: DenseVector[Double], k: Int): DenseMatrix[Double] = {
        val v = DenseMatrix.zeros[Double](dims, dims)

        v(ik(k),ik(k)) = -sin(theta(k))
    	v(ik(k),jk(k)) = cos(theta(k))
    	v(jk(k),ik(k)) = -cos(theta(k))
    	v(jk(k),jk(k)) = -sin(theta(k))

        return v
    }
}
