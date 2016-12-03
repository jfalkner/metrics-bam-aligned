package falkner.jayson.metrics.bam.aligned


import falkner.jayson.metrics.io.CSV
import org.specs2.mutable.Specification


class AlignedBamSpec extends Specification {

  "Aligned BAM" should {
    "Current version calculates without error" in {
      AlignedBam.version != null mustEqual true
    }
    "Support blank CSV generation" in {
      CSV(AlignedBam.blank).all != null mustEqual true
    }

  }
}
