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

import org.broadinstitute.gatk.engine.CommandLineGATK;

/**
 * Main class for GATKCT
 *
 * @author Daniel Gómez-Sánchez
 */
public class Main extends CommandLineGATK {

	public static void main(String[] argv) {
		printGATKCTheader();
		CommandLineGATK.main(argv);
	}

	public static void printGATKCTheader() {
		System.err.println("------------------------------------------------------------------------------------------");
		System.err.print("\t");
		System.err.println(ProjectProperties.getFormattedNameWithVersion());
		System.err.print("\t");
		System.err.println(ProjectProperties.getDescription());
		System.err.print("\t");
		System.err.println("Developed by " + ProjectProperties.getContact());
		System.err.println("------------------------------------------------------------------------------------------");
	}
}
