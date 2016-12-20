package falkner.jayson.metrics.bam.aligned

import java.nio.file.Paths

import falkner.jayson.metrics.io.JSON


object Main extends App {
    if (args.size != 2)
        println("Usage: java com.pacb.itg.metrics.bam.aligned.Main <file.bam> <output.json>")
    else
        JSON(Paths.get(args(1)), AlignedBam(Paths.get(args(0))))
}
