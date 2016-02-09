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
package org.magicdgs.gatkct;

import htsjdk.samtools.SAMException;
import htsjdk.tribble.TribbleException;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.exceptions.UserException;

/**
 * Main class for GATKCT
 *
 * @author Daniel Gómez-Sánchez
 */
public class Main extends CommandLineGATK {

	/**
	 * Main class copied from {@link org.broadinstitute.gatk.engine.CommandLineGATK} with addition of the custom tools
	 * at the begining and the end
	 *
	 * @param argv
	 */
	public static void main(String[] argv) {
		try {
			printGATKCTheader();
			CommandLineGATK instance = new CommandLineGATK();
			start(instance, argv);
			printGATKCend();
			System.exit(CommandLineProgram.result); // todo -- this is a painful hack
		} catch (UserException e) {
			exitSystemWithUserError(e);
		} catch (TribbleException e) {
			// We can generate Tribble Exceptions in weird places when e.g. VCF genotype fields are
			//   lazy loaded, so they aren't caught elsewhere and made into User Exceptions
			exitSystemWithUserError(e);
		} catch (SAMException e) {
			checkForMaskedUserErrors(e);
			exitSystemWithSamError(e);
		} catch (OutOfMemoryError e) {
			exitSystemWithUserError(new UserException.NotEnoughMemory());
		} catch (Throwable t) {
			checkForMaskedUserErrors(t);
			exitSystemWithError(t);
		}
	}

	/**
	 * Copied from {@link org.broadinstitute.gatk.engine.CommandLineGATK#checkForMaskedUserErrors(Throwable)}
	 *
	 * @param t the exception
	 */
	private static void checkForMaskedUserErrors(final Throwable t) {
		// masked out of memory error
		if (t instanceof OutOfMemoryError) {
			exitSystemWithUserError(new UserException.NotEnoughMemory());
		}
		// masked user error
		if (t instanceof UserException || t instanceof TribbleException) {
			exitSystemWithUserError(new UserException(t.getMessage()));
		}
		// no message means no masked error
		final String message = t.getMessage();
		if (message == null) {
			return;
		}
		// too many open files error
		if (message.contains("Too many open files")) {
			exitSystemWithUserError(new UserException.TooManyOpenFiles());
		}
		// malformed BAM looks like a SAM file
		if (message.contains(PICARD_TEXT_SAM_FILE_ERROR_1) || message.contains(PICARD_TEXT_SAM_FILE_ERROR_2)) {
			exitSystemWithSamError(t);
		}
		// can't close tribble index when writing
		if (message.contains("Unable to close index for")) {
			exitSystemWithUserError(new UserException(t.getCause() == null ? message : t.getCause().getMessage()));
		}
		// disk is full
		if (message.contains(NO_SPACE_LEFT_ON_DEVICE_ERROR) || message.contains(DISK_QUOTA_EXCEEDED_ERROR)) {
			exitSystemWithUserError(new UserException.NoSpaceOnDevice());
		}
		// masked error wrapped in another one
		if (t.getCause() != null) {
			checkForMaskedUserErrors(t.getCause());
		}
	}

	private static void printGATKCheaderSeparator() {
		System.err
			// .println("------------------------------------------------------------------------------------------");
			.println("==========================================================================================");
	}

	public static void printGATKCTheader() {
		printGATKCheaderSeparator();
		System.err.print("\t");
		System.err.print(ProjectProperties.getDescription());
		System.err.print(" (");
		System.err.print(ProjectProperties.getFormattedNameWithVersion());
		System.err.println(")");
		System.err.print("\t");
		System.err.println("Tools developed by " + ProjectProperties.getContact());
		printGATKCheaderSeparator();
	}

	public static void printGATKCend() {
		printGATKCheaderSeparator();
		System.err.print("\tThanks for using the tools implemented in ");
		System.err.println(ProjectProperties.getName());
		printGATKCheaderSeparator();
	}
}
