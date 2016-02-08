/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Daniel Gómez-Sánchez
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

package org.magicdgs.gatkct.tools;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import static org.broadinstitute.gatk.engine.SampleUtils.getUniqueSamplesFromRods;

/**
 * Mask genotypes for specific samples
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
 * @since 07-07-2015
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class MaskVariantsSampleSpecific extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {

	@ArgumentCollection
	protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

	@Output(doc="File to which variants should be written")
	protected VariantContextWriter writer = null;

	/**
	 * Any variant for the sample specified by -sn which overlaps entries from the provided mask rod will be call as missing.
	 * Note that there must be a 1-to-1 mapping between rod masks and sample names.
	 * If the user wants logic to be reversed, i.e. filter variants that do not overlap with provided mask, then argument
	 * -filterNotInMask can be used. Note that it is up to the user to adapt the name of the mask to make it clear that the
	 * reverse logic was used (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and
	 * -maskName=not_hapmap for the reverse masking).
	 */
	@Input(fullName="mask", shortName="mask", doc="Input sample specific ROD mask")
	public LinkedList<RodBinding<Feature>> mask;

	/**
	 * The sample_name that maps with the mask will be call as missing if overlaps with the provided mask.
	 * Note that there must be a 1-to-1 mapping between rod masks and sample names
	 */
	@Argument(fullName="sample_name", shortName="sn", doc="The sample to call as missing if 'mask' rod overlaps with a variant call")
	protected LinkedList<String> sampleNames;

	/**
	 * By default any variant falling in a mask will be filtered.
	 * If this argument is used, logic is reversed, and variants falling outside a given mask will be filtered.
	 * Use case is, for example, if we have an interval list or BED file with "good" sites.
	 * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
	 * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
	 */
	@Argument(fullName="filterNotInMask", shortName="filterNotInMask", doc="Mask sample genotypes NOT in given input mask.", required=false)
	protected boolean filterRecordsNotInMask = false;

	/**
	 * The default behavior of this tool is to remove only positions with missing genotypes for all the samples after masking.
	 * If this flag is set, all the variants will be written in the output regardless if the position is missing or not.
	 */
	@Argument(fullName="preserveAll", shortName="noRm", doc="No remove missing genotypes for all the samples")
	protected boolean preserveAll = false;

	/**
	 * By default, all the samples are called if they are not masked. If this argument is used, call as missing the genotypes for samples
	 * with less than the minimum coverage provided. If DP is not in a sample, it will be masked anyway.
	 */
	@Argument(fullName="minimum_coverage", shortName="minCov", doc="Minimum coverage for a sample to be called", required=false, minValue=0)
	protected int minCov = 0;

	// TODO: in development
	@Hidden
	@Argument(fullName="keepMaskedGT", shortName="keepGT", doc="Keep the masked GT in a FORMAT tag called "+MASKED_FORMAT_TAG, required=false)
	protected boolean keepMaskedGT=false;

	@Hidden
	@Argument(fullName="ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES", required=false, doc="Allow samples other than those in the VCF to be specified on the command line. These samples will be ignored.")
	private boolean allowNonOverlappingCommandLineSamples = false;

	//
	//	/**
	//	 * By default, all the samples are called if they are not masked. If this argument is used, call as missing the genotypes for samples
	//	 * with more than the maximum coverage provided.
	//	 */
	//	@Argument(fullName="maximum_coverage", shortName="maxCov", doc="Maximum coverage for a sample to be called", required=false, minValue=0)
	//	protected int maxCov = -1;
	//


	// Alleles for non-call genotypes
	private static final List<Allele> HAPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL);
	private static final List<Allele> DIPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
	// default value for previous genotype tag
	private static final String MASKED_FORMAT_TAG = "MGT";

	// HashMap to store the missing genotypes found in each sample. Maps sampleName -> number of missing
	private static ConcurrentHashMap<String, AtomicInteger> missingBySample;
	private static ConcurrentHashMap<String, RodBinding<Feature>> mapSampleMask;

	public void initialize() {
		// check if the masking files and names match
		if (mask.size() != sampleNames.size())
			throw new UserException.CommandLineException("--mask and --sample_name  must be a 1-to-1 mapping");
		// get the input names
		final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());
		Set<String> inputSamples = getUniqueSamplesFromRods(getToolkit(), inputNames);
		// initialize the missingBySample map
		missingBySample = new ConcurrentHashMap<>(sampleNames.size());
		mapSampleMask = new ConcurrentHashMap<>(sampleNames.size());
		// implementation for the ArrayList version
		//for(int i = 0; i < sampleNames.size(); i++) {
		//	String sample = sampleNames.get(i);
		//	RodBinding<Feature> masker = mask.get(i);
		while(!sampleNames.isEmpty()) {
			String sample = sampleNames.pop();
			RodBinding<Feature> masker = mask.pop();
			if(inputSamples.contains(sample)) {
				logger.info(String.format("Masking %s with positions %s %s file", sample, (filterRecordsNotInMask) ? "not in" : "in",  masker.getSource()));
				mapSampleMask.put(sample, masker);
				missingBySample.put(sample, new AtomicInteger(0));
			} else if(allowNonOverlappingCommandLineSamples) {
				logger.warn(String.format("Sample %s not found in the input file(s) and will be ignored", sample));
			} else {
				throw new UserException.BadInput(String.format("%s%s%n%n%s",
					"One sample entered on command line (through -sn) is not present in the VCF: ",
					sample,
					"To ignore this error, run with --allowNonOverlappingCommandLineSamples"));
			}
		}
		// genotypeFilterExps = VariantContextUtils.initializeMatchExps(genotypeFilterNames, genotypeFilterExpressions);

		VariantContextUtils.engine.get().setSilent(true);
		// setup the header fields
		Set<VCFHeaderLine> hInfo = new HashSet<>();
		hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));
		if(keepMaskedGT) {
			hInfo.add(GATKVCFHeaderLines.getFormatLine(MASKED_FORMAT_TAG));
		}
		writer.writeHeader(new VCFHeader(hInfo, inputSamples));
	}

	@Override
	public Integer reduceInit() {
		return 0;
	}

	@Override
	public Integer reduce(Integer value, Integer sum) {
		return value+sum;
	}

	@Override
	public Integer treeReduce(Integer lhs, Integer rhs) {
		return lhs + rhs;
	}

	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		if ( tracker == null )
			return 0;
		// get the variant context in the tracker
		Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
		// check if is null or empty
		if ( vcs == null || vcs.isEmpty()) {
			// logger.debug("This is returning a null boolean on"+context.getLocation());
			return 0;
		}
		int res = 0;
		// for each variant
		for(VariantContext vc: vcs) {
			// create a new variant
			VariantContextBuilder outputVariant = new VariantContextBuilder(vc);
			// obtain the new genotypes
			GenotypesContext outputGenotypes = callOverlapAsMissing(tracker, vc.getGenotypes(), context.getLocation());
			// check if it is null
			if(outputGenotypes != null) {
				// add the new genotypes
				outputVariant.genotypes(outputGenotypes);
				// write to the writer
				writer.add(outputVariant.make());
				res++;
			}
		}
		return res;
	}

	@Override
	public void onTraversalDone(Integer result) {
		logger.info(result + " records processed.");
		for(Map.Entry<String, AtomicInteger> a: missingBySample.entrySet()) {
			logger.info(String.format("%s genotypes called as missing for %s", a.getValue(), a.getKey()));
		}
	}

	/**
	 * Check if the sample specific ROD file overlaps with the position and masked if so. Return the new VariantContext
	 * with the masking or <code>null</code> if all are missing and preserveAll is not set
	 *
	 * @param tracker    the tracker for the variant
	 * @param original the variant context for all the samples
	 * @param loc	the location of the variant
	 * @return	the masked genotypes
	 */
	private GenotypesContext callOverlapAsMissing(RefMetaDataTracker tracker, GenotypesContext original, GenomeLoc loc) {
		int missingSamples = 0;
		// create the masked genotypes
		GenotypesContext maskedGenotypes = GenotypesContext.create();
		// for each of the original genotypes
		for(Genotype genotype: original) {
			// create the masked genotype as a builder
			GenotypeBuilder masked = new GenotypeBuilder(genotype);
			if(genotype.isCalled()) {
				if((genotype.getDP() < minCov)) {
					maskGenotypeBuilder(masked, genotype);
				} else {
					// get the ROD masker
					RodBinding<Feature> masker = mapSampleMask.get(genotype.getSampleName());
					// if there are a masker
					if (masker != null) {
						// check if the SNP is present in the corresponding mask
						boolean hasMask = (filterRecordsNotInMask) ? !tracker.hasValues(masker) : tracker.hasValues(masker);
						// change the masked genotype
						if (hasMask) {
							maskGenotypeBuilder(masked, genotype);
						}
					}
				}
			}
			Genotype newGenotype = masked.make();
			if(newGenotype.isNoCall()) {
				missingSamples++;
			}
			// add the masked genotype to the GenotypeContext to return
			maskedGenotypes.add(masked.make());
		}
		// if not all are preserved and all of them are missing, return null
		if( missingSamples >= original.size()) {
			if(preserveAll) {
				logger.warn("All missing genotypes preserved at "+ loc+".");
				return maskedGenotypes;
			} else {
				logger.debug("Missing genotypes at "+loc+".");
				return null;
			}
		}
		return maskedGenotypes;
	}

	private void maskGenotypeBuilder(GenotypeBuilder toMask, Genotype original) {
		missingBySample.get(original.getSampleName()).addAndGet(1);
		maskGenotypeBuilder(toMask, original.getPloidy());
		if(keepMaskedGT) {
			toMask.attribute(MASKED_FORMAT_TAG, original.getGenotypeString());
		}
	}


	/**
	 * Mask a genotype builder with the no_call alleles
	 *
	 * @param toMask	the genotype to mask
	 * @param ploidy	the ploidy to use
	 */
	private static void maskGenotypeBuilder(GenotypeBuilder toMask, int ploidy) {
		switch (ploidy) {
			case 1:
				toMask.alleles(HAPLOID_NO_CALL);
				break;
			case 2:
				toMask.alleles(DIPLOID_NO_CALL);
				break;
			default:
				toMask.alleles(Collections.nCopies(ploidy, Allele.NO_CALL));
				break;
		}
	}


	private static void addPreviousMaskedGenotype(GenotypeBuilder toMask, String originalGenotype) {
		toMask.attribute(MASKED_FORMAT_TAG, originalGenotype);
	}

}
