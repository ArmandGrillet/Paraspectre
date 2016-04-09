import breeze.linalg._
import breeze.numerics._
import breeze.stats._

package object serial {
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
        if (k >= distanceMatrix.cols - 1) {
            return max(distanceMatrix(*, ::)) // Maximum distance.
        } else {
            var localScale = DenseVector.zeros[Double](distanceMatrix.cols)
            var sortedVector = IndexedSeq(0.0)

            (0 until distanceMatrix.cols).map{col =>
                sortedVector = distanceMatrix(::, col).toArray.sorted
                localScale(col) = sortedVector(k) // Kth nearest distance.
            }

            return localScale
        }
    }

    def locallyScaledAffinityMatrix(distanceMatrix: DenseMatrix[Double], localScale: DenseVector[Double]): DenseMatrix[Double] = {
        var affinityMatrix = DenseMatrix.zeros[Double](distanceMatrix.rows, distanceMatrix.rows) // Distance matrix, size rows x rows.

        (0 until distanceMatrix.rows).map{ mainRow =>
            (mainRow + 1 until distanceMatrix.rows).map{ secondRow =>
                affinityMatrix(mainRow, secondRow) = -scala.math.pow(distanceMatrix(mainRow, secondRow), 2) // -d(si, sj)²
                affinityMatrix(mainRow, secondRow) = affinityMatrix(mainRow, secondRow) / (localScale(mainRow) * localScale(secondRow)) // -d(si, sj)² / lambi * lambj
                affinityMatrix(mainRow, secondRow) = scala.math.exp(affinityMatrix(mainRow, secondRow)) // exp(-d(si, sj)² / lambi * lambj)
                affinityMatrix(secondRow, mainRow) = affinityMatrix(mainRow, secondRow)
            }
        }

        return affinityMatrix
    }

    def logicalNot(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var result = DenseMatrix.zeros[Double](matrix.rows, matrix.cols)
        var row = 0
        while (row < matrix.rows) {
            var col = 0
            while (col < matrix.cols) {
                if (matrix(row, col) == 0.0) {
                    result(row, col) = 1.0
                }
                col += 1
            }
            row += 1
        }

        return result
    }

    def rotateEigenvectors(eigenvectors: DenseMatrix[Double]): (Double, DenseVector[Int], DenseMatrix[Double]) = {
        // Only works in two dimensions.
        val rows = eigenvectors.rows
        val cols = eigenvectors.cols

        // Get the number of angles
        val angles = (cols*(cols-1)/2).toInt
        val theta, newTheta = DenseVector.zeros[Double](angles)

        // We know that the method is 1

        // Build index mapping
        var ik = DenseVector.zeros[Int](angles)
        var jk = DenseVector.zeros[Int](angles)
        var k = 0
        var i = 0
        while (i < cols - 1) {
            var j = i + 1
            while (j < cols) {
                ik(k) = i
                jk(k) = j
                k += 1
                j += 1
            }
            i += 1
        }

        // Definitions
        val maxIterations = 200
        var iteration, angle = 0
        val stepSize = 1
        var old1cost, old2cost, cost, newCost = costForVectors(eigenvectors)

        while (iteration < maxIterations) {
            while (angle < angles) {
                val gradient = qualityGradient(eigenvectors, theta, ik, jk, angles, angle)
                newTheta(angle) = theta(angle) - stepSize * gradient
                val rotatedEigenvectors = rotateEigenvectorsWithGivenRotation(eigenvectors, newTheta, ik, jk, angles)
                newCost = costForVectors(rotatedEigenvectors)
                if (newCost < cost) {
                    theta(angle) = newTheta(angle)
                    cost = newCost
                } else {
                    newTheta(angle) = theta(angle)
                }
                angle += 1
            }
            if (iteration > 2 && abs(cost - old2cost) < 1e-3) { // Stopping criteria
                iteration = maxIterations
            } else {
                old2cost = old1cost
                old1cost = cost
            }
            iteration += 1
        }

        val rotatedEigenvectors = rotateEigenvectorsWithGivenRotation(eigenvectors, newTheta, ik, jk, angles)
        val clusts = clusters(rotatedEigenvectors)
        return (cost, clusts, rotatedEigenvectors)
    }

    def costForVectors(rotatedMatrix: DenseMatrix[Double]): Double = {
        // Take the square of all entries and find the max of each row
        val squaredMatrix = rotatedMatrix :* rotatedMatrix // :* = Hadamard product
        val maxRows = max(squaredMatrix(*, ::)) // We do not add a sqrt() as in the original code max_values[i] = p_X[ind]*p_X[ind];
        // Compute cost
        var cost = 0.0
        var row = 0
        while (row < rotatedMatrix.rows) {
            var col = 0
            while (col < rotatedMatrix.cols) {
                cost += squaredMatrix(row, col) / maxRows(row)
                col += 1
            }
            row += 1
        }

        return cost
    }

    def qualityGradient(eigenvectors: DenseMatrix[Double], theta: DenseVector[Double], ik: DenseVector[Int], jk: DenseVector[Int], angles: Int, index: Int): Double = {
        // In C++ it is a 1D array but arr(rows * row + col) = mat(row, col)
        var gradients = DenseMatrix.zeros[Double](eigenvectors.cols, eigenvectors.cols)
        gradients(ik(index), ik(index)) = -sin(theta(index))
        gradients(ik(index), jk(index)) = cos(theta(index))
        gradients(jk(index), ik(index)) = -cos(theta(index))
        gradients(jk(index), jk(index)) = -sin(theta(index))

        val u1 = uAB(theta, 0, index - 1, ik, jk, eigenvectors.cols)
        val u2 = uAB(theta, index + 1, angles - 1, ik, jk, eigenvectors.cols)

        val a = eigenvectors * gradients * u1 * u2

        // Rotate vectors according to current angles.
        val y = rotateEigenvectorsWithGivenRotation(eigenvectors, theta, ik, jk, angles)

        // Find the maximum of each row
        val squaredY = y :* y // :* = Hadamard product
        val maxEachRow = sqrt(max(squaredY(*, ::))) // Sqrt because in the original code max_values[i] = p_Y[ind];
        val argMaxEachRow = argmax(squaredY(*, ::)) // No sqrt because it is a position
        // Compute gradient
        var tmp1, tmp2, quality = 0.0
        var col = 0
        while (col < eigenvectors.cols) {
            var row = 0
            while (row < eigenvectors.rows) {
                tmp1 = a(row, col) * y(row, col) / scala.math.pow(maxEachRow(row), 2)
                tmp2 = a.toDenseVector(eigenvectors.rows * argMaxEachRow(row) + row) * squaredY(row, col) / scala.math.pow(maxEachRow(row), 3)
                quality += tmp1 - tmp2
                row += 1
            }
            col += 1
        }

        return 2 * quality / eigenvectors.rows / eigenvectors.cols
    }

    def uAB(theta: DenseVector[Double], a: Int, b: Int, ik: DenseVector[Int], jk: DenseVector[Int], cols: Int): DenseMatrix[Double] = {
        // Set uab to be an identity matrix
        var uab = DenseMatrix.eye[Double](cols)

        if (b < a) {
            return uab
        }

        var col = 0
        var k = a
        var tt, uIndex = 0.0
        while (k <= b) {
            tt = theta(k)
            while (col < cols) {
                uIndex = uab(ik(k), col) * cos(tt) - uab(jk(k), col) * sin(tt)
                uab(jk(k), col) = uab(ik(k), col) * sin(tt) * uab(jk(k), col) * cos(tt)
                uab(ik(k), col) = uIndex
                col += 1
            }
            k += 1
        }

        return uab
    }

    // Compute the Givens rotation
    // def givensRotation(n: Int, i: Int, j: Int, angle: Double): DenseMatrix[Double] = {
    //     return DenseMatrix.tabulate(n, n){
    //         case (row, col) =>
    //         if ((row == i && col == i) || (row == j && col == j)) {
    //             cos(angle)
    //         } else if (i > j && row == i && col == j) {
    //             sin(angle)
    //         } else if (row == j && col == i) {
    //             -sin(angle)
    //         } else if (row == col) {
    //             1.0
    //         } else {
    //             0.0
    //         }
    //     }
    // }

    def rotateEigenvectorsWithGivenRotation(eigenvectors: DenseMatrix[Double], theta: DenseVector[Double], ik: DenseVector[Int], jk: DenseVector[Int], angles: Int): DenseMatrix[Double] = {
        val g = uAB(theta, 0, angles - 1, ik, jk, eigenvectors.cols)
        val rotatedEigenvectors = eigenvectors * g
        return rotatedEigenvectors
    }

    def clusters(rotatedEigenvectors: DenseMatrix[Double]): DenseVector[Int] = {
        val squaredVectors = rotatedEigenvectors :* rotatedEigenvectors
        return argmax(squaredVectors(*, ::))
    }

    def positionInArray(value: Int, arr: Array[Int]): (Int, Int) = {
        var pos = 0
        for (pos <- (arr.length - 1) to 0 by -1) {
            if (value >= arr(pos)) {
                return (pos, arr(pos) - value)
            }
        }
        return (0, value)
    }
}
