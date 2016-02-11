/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Daniel Gómez-Sánchez
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.magicdgs.gatkct.tools.indels;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.*;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Hidden;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Identify regions with indels in BAM files
 *
 * <p>
 * [Functionality of this walker]
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 *
 * <h2>Examples</h2>
 * PRE-TAG
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 * PRE-TAG
 *
 * @author Daniel Gómez-Sánchez
 * @since 09-02-2016
 */
@DocumentedGATKFeature(groupName = HelpConstants.DOCS_CAT_DATA, extraDocs = {CommandLineGATK.class})
public class IdentifyIndelRegions extends LocusWalker<Interval, Integer>
	implements NanoSchedulable, TreeReducible<Integer> {

	/**
	 * The output intervals
	 */
	@Output(required = true)
	protected File out;

	/**
	 * The minimum count for indel identification. If the number of insertions or deletions of certain length at any
	 * position for all the inputs is lower than this number, the indels is not included in the output
	 */
	@Argument(fullName = "minimum-count", shortName = "m", doc = "Minimum count for an indel", required = false)
	int minCount = 1;

	/**
	 * Identified indel intervals will be extended this number of base-pairs in the sequence.
	 */
	@Argument(fullName = "indel-window", shortName = "w",
		doc = "Number of base-pairs to extend the window around the indel",
		required = false)
	int indelWin = 5;

	// this is the interval list that we will emit
	private IntervalList toEmit;

	private final static FormatUtil format = new FormatUtil();

	public void initialize() {
		super.initialize();
		// initialize the interval list with the header from the inputs
		toEmit = new IntervalList(getToolkit().getSAMFileHeader());
	}

	@Override
	public boolean includeReadsWithDeletionAtLoci() {
		// we should keep all the reads with an insertion
		return true;
	}

	/**
	 * Map the position to an interval if there is an indel
	 *
	 * @return the indel interval (without padding) if an indel is found; <code>null</code>
	 */
	@Override
	public Interval map(RefMetaDataTracker tracker, ReferenceContext reference, AlignmentContext alignment) {
		// get the pileup for this position
		ReadBackedPileup pileup = alignment.getBasePileup();
		// If there is no coverage, there is no indel region
		if (pileup.isEmpty()) {
			return null;
		}
		// initialize a new histogram
		final Histogram<Integer> delLegthCounts = new Histogram<>();
		int insCount = 0;
		// iterate over the pileup elements
		for (PileupElement element : pileup) {
			if (element.isDeletion()) {
				// if it is a deletion, add to the histogram of lengths
				delLegthCounts.increment(element.getCurrentCigarElement().getLength());
			} else if (element.isBeforeInsertion()) {
				// if it is an insertion
				insCount++;
			}
		}
		// if there is no indel, there are no region
		if (insCount == 0 && delLegthCounts.isEmpty()) {
			return null;
		}
		Interval indelInterval = new Interval(alignment.getContig(), (int) alignment.getPosition(),
			(int) alignment.getPosition());
		boolean emit = false;
		// iterate over each of the detected indels, to check if we should emit it and how much we did it
		for (Histogram<Integer>.Bin bin : delLegthCounts.values()) {
			if (bin.getValue() >= minCount) {
				emit = true;
				if (bin.getId()-1 > indelInterval.length()) {
					indelInterval = indelInterval.pad(0, bin.getId()-1);
				}
			}
		}
		// if it is not emited because of deletions, consider insertions
		if(!emit && insCount != 0) {
			if(insCount >= minCount) {
				emit = true;
				// create a 0-lenght interval
				indelInterval = indelInterval.pad(-1, 0);
			}
		}
		// emit the interval
		return (emit) ? indelInterval : null;
	}

	/**
	 * Add synchroniously the interval to the list of intervals
	 *
	 * @param interval the interval to add
	 */
	private synchronized void addToEmittedIntervals(Interval interval) {
		toEmit.add(interval);
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}

	@Override
	public Integer reduce(Interval interval, Integer value) {
		if (interval == null) {
			return value;
		}
		final Interval toAdd = interval.pad(indelWin, indelWin);
		// trying to add synchroniously
		synchronized (toEmit) {
			toEmit.add(toAdd);
		}
		// for each interval we will add one
		return value + 1;
	}

	@Override
	public Integer treeReduce(Integer lhs, Integer rhs) {
		return lhs + rhs;
	}

	public void onTraversalDone(Integer sum) {
		final NumberFormat fmt = new DecimalFormat("#,###");
		logger.info(String.format("Found %s positions with indels", fmt.format(sum)));
		toEmit = toEmit.uniqued();
		logger.info(String.format("Writting down the results in %s", out));
		long totalBpMasked = 0;
		final boolean intervalListFormat = FilenameUtils.getExtension(out.getName()).equals("interval_list");
		// write the header if it is an interval list format
		try (BufferedWriter bufferedWriter = IOUtil.openFileForBufferedWriting(out)) {
			if(intervalListFormat) {
				final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
				codec.encode(bufferedWriter, toEmit.getHeader());
			}
			for (Interval interval : toEmit) {
				totalBpMasked += interval.length();
				if (intervalListFormat) {
					writeIntervalListFormat(interval, bufferedWriter);
				} else {
					writeNotDefaultInterval(interval, bufferedWriter);
				}
			}
			bufferedWriter.flush();
			bufferedWriter.close();
		} catch (final IOException e) {
			throw new SAMException("Error writing out intervals to file: " + out.getAbsolutePath(), e);
		}
		logger.info(String.format("A total of %s intervals (%s bp) were identified", fmt.format(toEmit.size()), fmt.format(totalBpMasked)));
	}

	/**
	 * Write an interval in the output formatted as an interval list format
	 *
	 * @param interval the interval to write
	 * @param out      the writer where it should be write
	 *
	 * @throws IOException if there is a problem with the writer
	 */
	private void writeIntervalListFormat(Interval interval, BufferedWriter out) throws IOException {
		out.write(interval.getContig());
		out.write('\t');
		out.write(interval.getContig());
		out.write('\t');
		out.write(format.format(interval.getStart()));
		out.write('\t');
		out.write(format.format(interval.getEnd()));
		out.write('\t');
		out.write(interval.isPositiveStrand() ? '+' : '-');
		out.write('\t');
		if (interval.getName() != null) {
			out.write(interval.getName());
		} else {
			out.write(".");
		}
		out.newLine();
	}

	/**
	 * Write an interval in the output formatted as an non-default interval
	 *
	 * @param interval the interval to write
	 * @param out      the writer where it should be write
	 *
	 * @throws IOException if there is a problem with the writer
	 */
	private void writeNotDefaultInterval(Interval interval, BufferedWriter out) throws IOException {
		out.write(interval.getContig());
		out.write(':');
		out.write(format.format(interval.getStart()));
		if (interval.length() != 1) {
			out.write('-');
			out.write(format.format(interval.getEnd()));
		}
		out.newLine();
	}
}
