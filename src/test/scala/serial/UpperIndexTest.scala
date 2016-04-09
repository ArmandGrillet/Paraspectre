import serial._
import org.scalatest._

import breeze.linalg._
import breeze.numerics._
import breeze.stats._

class UpperIndexTest extends FlatSpec with Matchers {
    "The upper index " should "work with a 3x3 matrix and index = 2" in {
        val upperIdx = upperIndex(3, 2)
        val correctUpperIdx = (1, 2)
        upperIdx should be (correctUpperIdx)
    }
}
