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
package org.magicdgs.gatkct.tools.caller;

import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;

/**
 * Call variants in an sample that comes from ancient DNA
 * <p>
 * <p>
 * [Functionality of this walker]
 * </p>
 * <p>
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 * <p>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <p>
 * <h2>Examples</h2>
 * PRE-TAG
 * java
 * -jar GenomeAnalysisTK.jar
 * -T $WalkerName
 * PRE-TAG
 *
 * @author Daniel Gómez-Sánchez
 * @since 11-02-2016
 */
@DocumentedGATKFeature(groupName = HelpConstants.DOCS_CAT_VARDISC, extraDocs = {CommandLineGATK.class})
// TODO: create metrics class to be the reduce
public class AncientCaller extends LocusWalker<Integer, Long> {

    @Output(doc = "File to which variants should be written")
    protected VariantContextWriter writer = null;

    @Argument(fullName = "minimum_coverage", shortName = "minCov", doc = "Minimum coverage to call a variant",
            required = false)
    int minCov = 1;

    @Argument(fullName = "minimum_base_quality", shortName = "minBQ", doc = "Minimum base quality to use the bases",
            required = false)
    int minBQ = 1;

    @Argument(fullName = "minimum_mapping_quality", shortName = "minMQ",
            doc = "Minimum mapping quality for a read to be considered", required = false)
    int minMQ = 1;

    @Argument(fullName = "output_mode", shortName = "outMode",
            doc = "Specified which types of calls we should output", required = false)
    OutputOption outMode = OutputOption.CONFIDENT_VARIANTS;

    // TODO: add option to do not remove filtered bases

    // order of the bases in the base counts of
    private static final byte[] baseBytes = new byte[]{'A', 'C', 'T', 'G'};

    private String sampleName;

    private final static String LOW_COVERAGE_FILTER = "LowCov";
    private final static String POLYMORPHIC_FILTER = "Poly";

    public void initialize() {
        super.initialize();
        final GenomeAnalysisEngine toolkit = getToolkit();
        final Set<String> sampleNameSet = ReadUtils.getSAMFileSamples(toolkit.getSAMFileHeader());
        // TODO: make multi-sample?
        if (sampleNameSet.size() != 1) {
            String message = getClass().getSimpleName() + " only works with single ReadGroup BAM files";
            logger.error(message);
            throw new RuntimeException(message);
        } else {
            sampleName = sampleNameSet.iterator().next();
        }
        Set<VCFHeaderLine> headerSet = GATKVCFUtils.getHeaderFields(toolkit);
        // set info
        headerSet.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
        // set format
        headerSet.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        headerSet.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
        headerSet.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
        // TODO: I don't know if PASS should be added
        // headerSet.add(new VCFFilterHeaderLine(VCFConstants.PASSES_FILTERS_v4));
        // set filters
        if (!OutputOption.CONFIDENT_VARIANTS.equals(outMode)) {
            headerSet.add(new VCFFilterHeaderLine(LOW_COVERAGE_FILTER, "Coverage lower than " + minCov + " (user threshold)"));
            headerSet.add(new VCFFilterHeaderLine(POLYMORPHIC_FILTER, "Polymorphic site"));
        }
        writer.writeHeader(new VCFHeader(headerSet, new ArrayList<String>() {{
            add(sampleName);
        }}));
    }

    @Override
    public Integer map(RefMetaDataTracker refMetaDataTracker, ReferenceContext referenceContext,
                       AlignmentContext alignmentContext) {
        // Ns in the reference are not handle
        if (referenceContext.getBase() == 'N' || referenceContext.getBase() == 'n') {
            logger.debug("Found N at reference position" + referenceContext.getLocus());
            return 0; // we don't deal with the N ref base case
        }
        // Create the variant builder
        VariantContextBuilder variantBuilder = new VariantContextBuilder();
        variantBuilder.chr(referenceContext.getLocus().getContig());
        variantBuilder.loc(referenceContext.getLocus().getContig(),
                referenceContext.getLocus().getStart(),
                referenceContext.getLocus().getStop());

        // get the reference allele
        final Allele refAllele = Allele.create(referenceContext.getBase(), true);
        // get the pileup
        ReadBackedPileup pileup = alignmentContext.getBasePileup().getBaseAndMappingFilteredPileup(minBQ, minMQ);
        final int coverage = pileup.depthOfCoverage();
        variantBuilder.attribute(VCFConstants.DEPTH_KEY, coverage);
        if (coverage == 0) {
//			if (OutputOption.EMIT_ALL_SITES.equals(outMode)) {
//                variantBuilder.alleles(Collections.singleton(refAllele));
//                writer.add(variantBuilder.make());
//			}
            return 0;
        }
        Tuple<Genotype, Set<String>> callingResult = getGenotypeFromPileup(pileup, refAllele, minCov);
        if (callingResult == null) {
            return 0;
        }
        if(OutputOption.CONFIDENT_VARIANTS.equals(outMode) && !callingResult.b.isEmpty()) {
            return 0;
        }
        // TODO: this should be changed for the multi-sample implementation
        // creating the allele set for the variant builder
        HashSet<Allele> alleleSet = new HashSet<>(callingResult.a.getAlleles());
        alleleSet.add(refAllele);
        variantBuilder.alleles(alleleSet);
        variantBuilder.genotypes(callingResult.a);
        // TODO: check if it output PASS
        variantBuilder.filters(callingResult.b);
        writer.add(variantBuilder.make(true));
        return 1;
    }

    /**
     * Get the genotype for a concrete sample from the pileup
     *
     * @param pileup    the pileup for the sample
     * @param refAllele the reference allele in this position
     * @param minCov    the minimum coverage to do not filter
     * @return the genotype and the filters
     * @throws IllegalArgumentException if the pileup contains more than one sample
     */
    private static Tuple<Genotype, Set<String>> getGenotypeFromPileup(final ReadBackedPileup pileup, final Allele refAllele, int minCov) {
        final Collection<String> samples = pileup.getSamples();
        final HashSet<String> filter = new LinkedHashSet<>();
        if (samples.size() != 1) {
            throw new IllegalArgumentException("Pileup should be sample specific");
        }
        String sampleName = samples.iterator().next();
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final int coverage = pileup.depthOfCoverage();
        genotypeBuilder.DP(coverage);
        if (coverage == 0) {
            // return a genotype without anything
            return new Tuple<>(genotypeBuilder.make(), filter);
        } else if (coverage < minCov) {
            filter.add(LOW_COVERAGE_FILTER);
        }
        // get the base count in the order A, C, G, T
        int[] baseCounts = pileup.getBaseCounts();
        // get the reference base
        final HashMap<Allele, Integer> alleleCounts = new HashMap<>(2);
        for (int i = 0; i < baseCounts.length; i++) {
            if (baseCounts[i] != 0) {
                // add the alleles
                alleleCounts.put(Allele.create(baseBytes[i], baseBytes[i] == refAllele.getBases()[0]),
                        baseCounts[i]);
            }
        }
        final ArrayList<Allele> alleles = new ArrayList<>(alleleCounts.keySet());
        switch (alleles.size()) {
            case 0: // this is not possible, because the depth is checked previously
                logger.error("Unreacheable code");
                throw new RuntimeException("Unreacheable code");
            case 1: // monomorphic site
                // duplicate the allele
                alleles.add(alleles.get(0));
                break;
            case 2: // polymorphic sites
                genotypeBuilder.alleles(alleles);
                filter.add(POLYMORPHIC_FILTER);
            default: // TODO: warning because it should be diploid
                return null;
        }
        // setting alleles
        genotypeBuilder.alleles(alleles);
        final Integer refCount = alleleCounts.getOrDefault(refAllele, 0);
        alleleCounts.remove(refAllele);
        final Integer alt = (alleleCounts.isEmpty()) ? 0 : alleleCounts.values().stream().findFirst().get();
        // because only biallelic
        genotypeBuilder.AD(new int[]{refCount, alt});
        // creating the allele set for the variant builder
        return new Tuple<>(genotypeBuilder.make(), filter);
    }

    @Override
    public Long reduceInit() {
        return 0L;
    }

    @Override
    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }
}
