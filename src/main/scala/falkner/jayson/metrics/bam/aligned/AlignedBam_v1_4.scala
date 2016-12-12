package falkner.jayson.metrics.bam.aligned

import java.nio.file.{Files, Path}
import java.text.DecimalFormat

import falkner.jayson.metrics.Distribution._
import falkner.jayson.metrics._
import falkner.jayson.metrics.io.CSV
import htsjdk.samtools.{SAMFileHeader, SAMRecord, SamReaderFactory, ValidationStringency}

import scala.collection.JavaConverters._
import scala.concurrent.{Await, Future}
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.Duration
import scala.util.{Failure, Success, Try}

object AlignedBam_v1_4 {
  case class Read(name: String, len: Int, lenCovered: Int, refLen: Int, mappingQuality: Int, refAccClip: Float, refAccNoClip: Float, readAccClip: Float, readAccNoClip: Float, cigar: Cigar)
  case class Cigar(M: Int, I: Int, D: Int, N: Int, S: Int, H: Int, P: Int, EQ:Int, X: Int)

  lazy val sixFractionDigits = new DecimalFormat("#.######")

  class ReadMetric(r: Read) extends Metrics {
    override val namespace: String = "Read"
    override val version: String = "_"
    override val values: List[Metric] = List(
      Str("Name", r.name),
      Num("Read Length", r.len),
      Num("Ref Length", r.refLen),
      Num("Mapping Quality", r.mappingQuality, sixFractionDigits),
      Num("Ref Accuracy (Clip)", r.refAccClip, sixFractionDigits),
      Num("Ref Accuracy (No Clip)", r.refAccNoClip, sixFractionDigits),
      Num("Read Accuracy (Clip)", r.readAccClip, sixFractionDigits),
      Num("Read Accuracy (No Clip)", r.readAccNoClip, sixFractionDigits),
      Num("M", r.cigar.M),
      Num("D", r.cigar.D),
      Num("I", r.cigar.I),
      Num("N", r.cigar.N),
      Num("S", r.cigar.S),
      Num("H", r.cigar.H),
      Num("P", r.cigar.P),
      Num("=", r.cigar.EQ),
      Num("X", r.cigar.X)
    )
  }

  def apply(p: Path): AlignedBam_v1_4 = new AlignedBam_v1_4(p)

  val version: String = "1.4"

  def exportReads(reads: Seq[Read], p: Path): Unit = {
    // TODO: auto-delete on error?
    val bw = Files.newBufferedWriter(p)
    Seq(reads.head).map(r => bw.write(CSV(new ReadMetric(r)).all + "\n"))
    reads.tail.foreach(r => bw.write(CSV(new ReadMetric(r)).values + "\n"))
    bw.flush()
    bw.close()
  }

  case class Chunk(size: Long,
                   totalReadLength: Long,
                   //totalReadLengthCovered: Long,
                   // cigar op tallies
                   totalMatchOrMismatch: Long,
                   totalIns: Long,
                   totalDel: Long,
                   totalSkip: Long,
                   totalSoftClip: Long,
                   totalHardClip: Long,
                   totalPadding: Long,
                   totalMatch: Long,
                   totalMismatch: Long,
                   // distributions
                   readLength: Discrete,
                   mappingQuality: Discrete,
                   refAccuracyWithClip: Continuous,
                   refAccuracyNoClip: Continuous,
                   readAccuracyWithClip: Continuous,
                   readAccuracyNoClip: Continuous,
                  // cigar op frequency dists
                   matchOrMismatchFrequency: Continuous,
                   insFrequency: Continuous,
                   delFrequency: Continuous,
                   skipFrequency: Continuous,
                   softClipFrequency: Continuous,
                   hardClipFrequency: Continuous,
                   paddingFrequency: Continuous,
                   matchFrequency: Continuous,
                   mismatchFrequency: Continuous)

  val (freqMin, freqMax) = (0.0f, 1.0f)
}

/**
  * Companion class that calculates metrics via single-pass through a BAM file
  *
  * See README.md for details about each metric.
  */
class AlignedBam_v1_4(p: Path, nBins: Int = 30) extends Metrics {
  import AlignedBam_v1_4._
  override val namespace: String = "BAM"
  override val version: String = s"${AlignedBam.version}~${AlignedBam_v1_4.version}"
  override val values: List[Metric] = List(
    Str("Code Version", AlignedBam.version),
    Str("Spec Version", AlignedBam_v1_4.version),
    Num("Reads", chunks.map(_.size).sum),
    Num("Total Read Length", chunks.map(_.totalReadLength).sum),
    // totals for each cigar element
    Num("Total Match or Mismatch", chunks.map(_.totalMatchOrMismatch).sum),
    Num("Total Ins", chunks.map(_.totalIns).sum),
    Num("Total Del", chunks.map(_.totalDel).sum),
    Num("Total Skip", chunks.map(_.totalSkip).sum),
    Num("Total Soft Clip", chunks.map(_.totalSoftClip).sum),
    Num("Total Hard Clip", chunks.map(_.totalHardClip).sum),
    Num("Total Padding", chunks.map(_.totalHardClip).sum),
    Num("Total Match", chunks.map(_.totalMatch).sum),
    Num("Total Mismatch", chunks.map(_.totalMismatch).sum),
    //Num("Total Read Length (Covered)", chunks.map(_.totalReadLengthCovered).sum),
    Dist("Read Length", mergeDiscrete(chunks.map(_.readLength))),
    Dist("Mapping Quality", mergeDiscrete(chunks.map(_.mappingQuality))),
    // accuracies
    DistCon("Accuracy (Clip)", mergeContinuous(chunks.map(_.refAccuracyWithClip), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Accuracy (No Clip)", mergeContinuous(chunks.map(_.refAccuracyNoClip), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Accuracy (Read-based w/Clip)", mergeContinuous(chunks.map(_.readAccuracyWithClip), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Accuracy (Read-based No Clip)", mergeContinuous(chunks.map(_.readAccuracyNoClip), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    // mutation frequencies
    DistCon("Match or Mismatch Freq.", mergeContinuous(chunks.map(_.matchOrMismatchFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Del Freq.", mergeContinuous(chunks.map(_.delFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Ins Freq.", mergeContinuous(chunks.map(_.insFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Skip Freq.", mergeContinuous(chunks.map(_.skipFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Soft Clip Freq.", mergeContinuous(chunks.map(_.softClipFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Hard Clip Freq.", mergeContinuous(chunks.map(_.hardClipFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Padding Freq.", mergeContinuous(chunks.map(_.paddingFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Match Freq.", mergeContinuous(chunks.map(_.matchFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits),
    DistCon("Mismatch Freq.", mergeContinuous(chunks.map(_.mismatchFrequency), forceMin = Some(freqMin), forceMax = Some(freqMax)), sixFractionDigits)
  )

  def handleReads(buf: Seq[Read]): Chunk = Chunk(
    buf.size,
    buf.map(_.len).sum, // read length
    // mutation tallies
    buf.map(_.cigar.M.toLong).sum,
    buf.map(_.cigar.I.toLong).sum,
    buf.map(_.cigar.D.toLong).sum,
    buf.map(_.cigar.N.toLong).sum,
    buf.map(_.cigar.S.toLong).sum,
    buf.map(_.cigar.H.toLong).sum,
    buf.map(_.cigar.P.toLong).sum,
    buf.map(_.cigar.EQ.toLong).sum,
    buf.map(_.cigar.X.toLong).sum,
    // distributions
    calcDiscrete(buf.map(_.len)),
    calcDiscrete(buf.map(_.mappingQuality)), // mappingQualityDist,
    // accuracies
    calcContinuous(buf.map(_.refAccClip), forceMin = Some(freqMin), forceMax = Some(freqMax)),
    calcContinuous(buf.map(_.refAccNoClip), forceMin = Some(freqMin), forceMax = Some(freqMax)),
    calcContinuous(buf.map(_.readAccClip), forceMin = Some(freqMin), forceMax = Some(freqMax)),
    calcContinuous(buf.map(_.readAccNoClip), forceMin = Some(freqMin), forceMax = Some(freqMax)),
    // mutation frequencies
    calcContinuous(buf.map(v => v.cigar.M.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // match or mismatch frequency
    calcContinuous(buf.map(v => v.cigar.I.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // ins frequency
    calcContinuous(buf.map(v => v.cigar.D.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // del frequency
    calcContinuous(buf.map(v => v.cigar.N.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // skip frequency
    calcContinuous(buf.map(v => v.cigar.S.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // soft clip Frequency
    calcContinuous(buf.map(v => v.cigar.H.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // hard clip Frequency
    calcContinuous(buf.map(v => v.cigar.P.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // padding Frequency
    calcContinuous(buf.map(v => v.cigar.EQ.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax)), // match frequency
    calcContinuous(buf.map(v => v.cigar.X.toFloat / v.len), forceMin = Some(freqMin), forceMax = Some(freqMax))) // mismatch frequency)

  def makeRead(r: SAMRecord): Read = {
    // calculate all per-read metrics in parallel
    val readLen = r.getReadLength
    val refLen = r.getCigar.getReferenceLength
    val c = r.getCigar.getCigarElements.asScala.toList.groupBy(_.getOperator.toString).toList.map { case (k, v) => (k, v.map(_.getLength).sum) }.toMap
    Read(
      r.getReadName,
      readLen,
      r.getCigar.getReadLength,
      refLen,
      r.getMappingQuality,
      refAccuracyWithClip(refLen, c),
      refAccuracyNoClip(refLen, c),
      readAccuracyWithClip(readLen, c),
      readAccuracyNoClip(readLen, c),
      Cigar(
        c.getOrElse("M" , 0),
        c.getOrElse("I" , 0),
        c.getOrElse("D" , 0),
        c.getOrElse("N" , 0),
        c.getOrElse("S" , 0),
        c.getOrElse("H" , 0),
        c.getOrElse("P" , 0),
        c.getOrElse("=" , 0),
        c.getOrElse("X" , 0)
      )
    )
  }

  lazy val chunkSize = 10000

  lazy val (header, chunks): (SAMFileHeader, List[Chunk]) = Try {
    val factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
    val bam = factory.open(p)

    (bam.getFileHeader, bam.iterator.asScala.grouped(chunkSize).map(g =>
      handleReads(g.map(r => Future(makeRead(r))).map(fr => Await.result(fr, Duration.Inf)))).toList)
  } match {
    case Success(s) =>
      println("Done processing.")
      println(s"Made ${s._2.size} chunks")
      s
    case Failure(t) if p == null => (null, null) // support AlignedPacBioBam.blank
    case Failure(t) => throw t
  }

  def accuracy(len: Int, keys: Seq[String], cigarMap: Map[String, Int]): Float =
    return 1f - (keys.flatMap(cigarMap.get).sum.toFloat / len)

  def accuracyWithClip(l: Int, c: Map[String, Int]): Float = accuracy(l, Seq("X", "I", "D", "N", "S", "H"), c)

  def accuracyNoClip(l: Int, c: Map[String, Int]): Float = accuracy(l, Seq("X", "I", "D"), c)

  def refAccuracyWithClip(refLen: Int, cigarMap: Map[String, Int]): Float = accuracyWithClip(refLen, cigarMap)

  def refAccuracyNoClip(refLen: Int, cigarMap: Map[String, Int]): Float = accuracyNoClip(refLen, cigarMap)

  def readAccuracyWithClip(readLen: Int, cigarMap: Map[String, Int]): Float = accuracyWithClip(readLen, cigarMap)

  def readAccuracyNoClip(readLen: Int, cigarMap: Map[String, Int]): Float = accuracyNoClip(readLen, cigarMap)

  // could dump stats for base call transition frequency per (1 prior base, 2 prior bases) ins and del -- i.e. does enzyme mess up during transition

  // dump stats for incorrect base call per sequence location -- i.e. any bias toward beginning, end, etc.

  // dump stats for base quality vs incorrect call -- i.e. do we need to better threshold on quality or otherwise split reads?
}