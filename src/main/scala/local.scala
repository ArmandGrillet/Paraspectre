import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import java.io.File

object Local {
    def main(args: Array[String]) = {
        // Choose the dataset to cluster.
        val pathToMatrix = getClass.getResource("/0.csv").getPath()
        val matrixFile = new File(pathToMatrix)

        // Create a DenseMatrix from the CSV
        var matrix = breeze.linalg.csvread(matrixFile)
        // println(matrix)

        // Centralizing and scala the data
        val meancols = mean(matrix(::, *))
        // Waiting for fix scalanlp/breeze#450
        matrix = (matrix.t(::,*) - meancols.t).t
        matrix /= max(abs(matrix))

        println(matrix)
    }

    def rowrepmat(matrix: breeze.linalg.DenseMatrix[Double], rep: Int): breeze.linalg.DenseMatrix[Double] = {
        var X = matrix.copy
        for(a <- 0 to rep){
            X = DenseMatrix.vertcat(X, matrix)
        }
        return X
    }

    // Compute the mean of each column
    def meancols(matrix: breeze.linalg.DenseMatrix[Double]): Array[String] = {
        var mean = new Array[String](matrix.cols)
        return mean
    }

}
