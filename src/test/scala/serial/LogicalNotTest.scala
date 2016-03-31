import serial._
import org.scalatest._

import breeze.linalg._
import breeze.numerics._
import breeze.stats._

class LogicalNotTest extends FlatSpec with Matchers {
    "The logical not " should "work with a regular matrix" in {
        val matrix = DenseMatrix((-253.0, 0.0), (12.0, -2.0), (2.0, 0.0))
        val logicNot = logicalNot(matrix)
        val correctLogicNot = DenseMatrix((0.0, 1.0), (0.0, 0.0), (0.0, 1.0))
        logicNot should be (correctLogicNot)
    }

    "The logical not " should "work with a big zero matrix" in {
        val matrix = DenseMatrix.zeros[Double](999, 999)
        val logicNot = logicalNot(matrix)
        val correctLogicNot = DenseMatrix.ones[Double](999, 999)
        logicNot should be (correctLogicNot)
    }

    "The logical not " should "work with a big one matrix" in {
        val matrix = DenseMatrix.ones[Double](999, 999)
        val logicNot = logicalNot(matrix)
        val correctLogicNot = DenseMatrix.zeros[Double](999, 999)
        logicNot should be (correctLogicNot)
    }
}
