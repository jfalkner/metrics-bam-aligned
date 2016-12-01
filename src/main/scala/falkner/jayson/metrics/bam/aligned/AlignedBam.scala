package falkner.jayson.metrics.bam.aligned

import java.nio.file.Path

object AlignedBam {

  def apply(p: Path): AlignedBam_v1_4 = new AlignedBam_v1_4(p)

  val version: String = "0.0.1"
}