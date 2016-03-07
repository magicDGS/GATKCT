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
package org.magicdgs.gatkct.tools.misc;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * For each sample, assign randomly each allele to each chromosome, and obtain one sample per chromosome.
 * This walker only works with diploid individuals
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
 * @since 07-07-2015
 */
public class GetRandomChromosomes extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc = "File to which variants should be written")
    protected VariantContextWriter writer = null;

    @Argument(fullName = "seed", doc = "Seed for randomize the chromosomes", required = false)
    public long seed = 1;

    /**
     * This argument can be specified multiple times in order to provide multiple sample names.
     */
    @Argument(fullName = "sample_name", shortName = "sn", doc = "Include this sample", required = false)
    public Set<String> sampleNames = new HashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. This
     * argument can be specified multiple times in order to use multiple different matching patterns.
     */
    @Argument(fullName = "sample_expressions", shortName = "se", doc = "Regular expression to select multiple samples", required = false)
    public Set<String> sampleExpressions;

    /**
     * Sample names should be in a plain text file listing one sample name per line. This argument can be specified multiple times in order to provide
     * multiple sample list files.
     */
    @Input(fullName = "sample_file", shortName = "sf", doc = "File containing a list of samples to include", required = false)
    public Set<File> sampleFiles;

    /**
     * Note that sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be
     * excluded. This argument can be specified multiple times in order to provide multiple sample names.
     */
    @Argument(fullName = "exclude_sample_name", shortName = "xl_sn", doc = "Exclude this sample", required = false)
    public Set<String> XLsampleNames = new HashSet<>(0);

    /**
     * Sample names should be in a plain text file listing one sample name per line. Note that sample exclusion takes precedence over inclusion, so that
     * if a sample is in both lists it will be excluded. This argument can be specified multiple times in order to
     * provide multiple sample list files.
     */
    @Input(fullName = "exclude_sample_file", shortName = "xl_sf", doc = "List of samples to exclude", required = false)
    public Set<File> XLsampleFiles = new HashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. Note that sample exclusion takes precedence
     * over inclusion, so that if a sample is in both lists it will be excluded. This  argument can be specified multiple times in order to use multiple
     * different matching patterns.
     */
    @Input(fullName = "exclude_sample_expressions", shortName = "xl_se", doc = "List of sample expressions to exclude", required = false)
    public Set<String> XLsampleExpressions = new HashSet<>(0);

    // the samples to include
    private TreeSet<String> samples = new TreeSet<>();

    private Random randomizer;

    public void initialize() {
        // Get list of samples to include in the output
        // first get the sample for the RODs
        // Get list of samples to include in the output
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<>(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        // get the samples for the files
        Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
        // get the samples for the expression
        Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

        // first, check overlap between requested and present samples
        Set<String> commandLineUniqueSamples = new HashSet<>(samplesFromFile.size() + samplesFromExpressions.size() + sampleNames.size());
        commandLineUniqueSamples.addAll(samplesFromFile);
        commandLineUniqueSamples.addAll(samplesFromExpressions);
        commandLineUniqueSamples.addAll(sampleNames);
        commandLineUniqueSamples.removeAll(vcfSamples);

        // second, add the requested samples
        samples.addAll(sampleNames);
        samples.addAll(samplesFromExpressions);
        samples.addAll(samplesFromFile);

        logger.debug(Utils.join(",", commandLineUniqueSamples));
        boolean noSamplesSpecified = false;
        // if none were requested, we want all of them
        if (samples.isEmpty()) {
            samples.addAll(vcfSamples);
            noSamplesSpecified = true;
        }
        // now, exclude any requested samples
        final Collection<String> XLsamplesFromFile = SampleUtils.getSamplesFromFiles(XLsampleFiles);
        final Collection<String> XLsamplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, XLsampleExpressions);
        samples.removeAll(XLsamplesFromFile);
        samples.removeAll(XLsampleNames);
        samples.removeAll(XLsamplesFromExpressions);
        noSamplesSpecified = noSamplesSpecified && XLsampleNames.isEmpty() && XLsamplesFromFile.isEmpty() &&
                XLsamplesFromExpressions.isEmpty();

        if (samples.isEmpty() && !noSamplesSpecified) {
            throw new UserException("All samples requested to be included were also requested to be excluded.");
        }
        // log the samples included and generate the samples for the writer
        HashSet<String> outputSamples = new HashSet<>(samples.size() * 2);
        for (String sample : samples) {
            if (!noSamplesSpecified) {
                logger.info("Including sample '" + sample + "'");
            }
            outputSamples.add(getSampleChromosomeName(sample, 1));
            outputSamples.add(getSampleChromosomeName(sample, 1));
        }
        // Initialize VCF header
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFHeaderLine("source", this.getClass().getSimpleName()));
        // write the header with the output samples
        writer.writeHeader(new VCFHeader(headerLines, outputSamples));
        // set the randomizer
        randomizer = new Random(seed);
    }

    /**
     * Get the chromosome name for this sample
     *
     * @param sample the sample name
     * @param chr    the number of the chromosome
     * @return the name as sample_n
     */
    private static String getSampleChromosomeName(String sample, int chr) {
        return String.format("%s_%s", sample, chr);
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return 0;
        // get the variant context in the tracker
        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
        // check if is null or empty
        if (vcs == null || vcs.isEmpty()) {
            return 0;
        }
        int res = 0;
        // for each variant
        for (VariantContext vc : vcs) {
            // create a new variant and add to the writer
            writer.add(new VariantContextBuilder(vc)
                    .genotypes(createRandomChromosomes(vc.getGenotypes(samples)))
                    .make()
            );
            res++;
        }
        return res;
    }

    /**
     * Create random chromosomes for the genotypes
     *
     * @param genotypes the genotypes
     * @return the duplicated genotypes
     */
    private GenotypesContext createRandomChromosomes(final GenotypesContext genotypes) {
        GenotypesContext toReturn = GenotypesContext.create(genotypes.size() * 2);
        for (Genotype geno : genotypes) {
            // check if the sample is diploid
            if (geno.getPloidy() != 2) {
                throw new UserException(this.getClass().getSimpleName() + " found a non-diploid sample.");
            }
            if (geno.isHom()) {
                // easy to duplicate if it is homozygous
                final String sampleName = geno.getSampleName();
                final GenotypeBuilder builder = new GenotypeBuilder(geno);
                toReturn.add(builder.name(getSampleChromosomeName(sampleName, 1)).make());
                toReturn.add(builder.name(getSampleChromosomeName(sampleName, 2)).make());
            } else {
                // get a random sample
                final int random = randomizer.nextInt(1);
                toReturn.add(createGenotypeWithAllele(geno, 1, random));
                toReturn.add(createGenotypeWithAllele(geno, 2, Math.abs(random - 1)));
            }

        }
        return toReturn;
    }

    /**
     * Get the genotype for the chromosome specified
     *
     * @param genotype    the genotype for the sample
     * @param chr         the number of the chromosome
     * @param alleleIndex the index of the allele to use
     * @return homozygous genotype for the chromosome with the allele in the specified index
     */
    private Genotype createGenotypeWithAllele(Genotype genotype, int chr, int alleleIndex) {
        return new GenotypeBuilder(genotype)
                .name(getSampleChromosomeName(genotype.getSampleName(), chr))
                .alleles(Collections.nCopies(2, genotype.getAllele(alleleIndex)))
                .make();
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }


    @Override
    public void onTraversalDone(Integer result) {
        logger.info(result + " variants processed.");
    }
}
