package serial

import breeze._
import java.io.File

object Serial {
    // The CLI tool
    def main(args: Array[String]) {
        var dataset = getClass.getResource("/0.csv").getPath()
        var min = 2
        var max = 6
        var debug = false
        var result = false

        args.sliding(2, 1).toList.collect {
            case Array("--data", argData: String) => dataset = argData
            case Array("--sample", argSample: String) => dataset = getClass.getResource("/" + argSample + ".csv").getPath()
            case Array("--min", argMin: String) => min = argMin.toInt
            case Array("--max", argMax: String) => max = argMax.toInt
            case Array("--debug", argDebug: String) => debug = argDebug.toBoolean
            case Array("--result", argResult: String) => result = argResult.toBoolean
        }

        if (dataset.takeRight(4) != ".csv") {
            println("Unusable dataset, needs a .csv file.")
            return
        }

        if (min > max) {
            println(min + " > " + max)
            return
        }

        // Load the dataset to cluster.
        val datasetFile = new File(dataset)

        // Create a DenseMatrix from the CSV.
        val originalMatrix = breeze.linalg.csvread(datasetFile)

        val algorithm = new Algorithm(min, max, debug)
        val clusters = algorithm.cluster(originalMatrix)
        if (result) {
            println(clusters)
        }
    }
}
