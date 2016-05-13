package serial

import breeze._
import java.awt.image._
import java.io.File

object Serial {
    // The CLI tool
    def main(args: Array[String]) {
        var dataPath = ""
        var min = 2
        var max = 6
        var debug = false
        var result = false

        args.sliding(2, 1).toList.collect {
            case Array("--data", argData: String) => dataPath = argData
            case Array("--sample", argSample: String) => dataPath = getClass.getResource("/" + argSample + ".csv").getPath()
            case Array("--min", argMin: String) => min = argMin.toInt
            case Array("--max", argMax: String) => max = argMax.toInt
            case Array("--debug", argDebug: String) => debug = argDebug.toBoolean
            case Array("--result", argResult: String) => result = argResult.toBoolean
        }

        if (min > max) {
            println(min + " > " + max)
            return
        }

        if (dataPath.takeRight(4) != ".csv") {
            println("Unusable dataset, needs a .csv file.")
            return
        }
        // Load the dataset to cluster.
        val dataset = new File(dataPath)

        // Create a DenseMatrix from the CSV.
        val matrix = breeze.linalg.csvread(dataset)

        val algorithm = new Algorithm(min, max, debug)
        val clusters = algorithm.cluster(matrix)
        if (result) {
            println(clusters)
        }
    }
}
