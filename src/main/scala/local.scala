import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import java.io.File

object Local {
    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/0.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV.
        var matrix = breeze.linalg.csvread(matrixFile)
        // println(matrix)

        // Centralizing and scale the data.
        val meancols = mean(matrix(::, *))
        // Waiting for fix scalanlp/breeze#450
        matrix = (matrix.t(::, *) - meancols.t).t
        matrix /= max(abs(matrix))

        // Build affinity matrix.
        val distances = distanceMatrix(matrix)

        // println(distances)
    }

    def distanceMatrix(matrix: DenseMatrix[Double]): DenseMatrix[Double] = {
        var distanceMatrix = DenseMatrix.zeros[Double](matrix.rows,matrix.rows)
        var distanceVector = DenseVector(0.0).t
        var distance = 0.0
        (0 until matrix.rows).map{mainRow =>
            (mainRow + 1 until matrix.rows).map{secondRow =>
                distanceVector = matrix(mainRow, ::) - matrix(secondRow,::) // Xi - Xj
                distanceVector *= distanceVector // (Xi - Xj)^^2
                distance = sqrt(sum(distanceVector)) // √(Xi - Xj)^^2 + (Yi - Yj)^^2 + ...
                distanceMatrix(mainRow, secondRow) = distance
                distanceMatrix(secondRow, mainRow) = distance
            }
        }
        println(distanceMatrix)
        return distanceMatrix
    }
}
