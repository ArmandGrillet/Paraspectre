import serial._
import org.scalatest._

import breeze.linalg._
import breeze.numerics._
import breeze.stats._

class CostTest extends FlatSpec with Matchers {
    "The cost " should "work with a perfect 2x2 matrix" in {
        val matrix = DenseMatrix((1.0, 0.0), (0.0, 1.0))
        val computedCost = cost(matrix)
        val correctCost = 2.0
        computedCost should be (correctCost)
    }
}
