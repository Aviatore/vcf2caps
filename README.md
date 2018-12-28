# VCF2CAPS

VCF2CAPS is a VCF file analysis tool for CAPS marker development. The software facilitates the conversion of a large number of single nucleotide polymorphisms (SNPs), multiple nucleotide polymorphisms (MNPs) and insertion/deletion (InDel) polymorphisms detected by SNP calling tools, e.g. SAMtools, Platypus or FreeBayes, into CAPS markers.

## Getting started

The software is provided as a Perl script 'vcf2caps.pl' which needs Perl environment to be installed to run properly. However, if you are working on Windows OS you do not need to install anything. Simply run the executable file VCF2CAPS.exe.
The detailed description of how to use this software is available in the manual file.

### Prerequisites

- Perl v.5.008 or later – the Perl language interpreter.
- Perl modules used by vcf2caps.pl:
  - Tk
  - Tk::Notebook
  - Tk::FBox
  - Tk::Animation
  - Tk::Listbox
  - Tk::Pane
  - Tk::ProgressBar
  - Tk::PNG
  - LWP::Simple
  - Digest::MD5
  - threads
  - threads::shared
  - Encode
  - utf8

### Install dependencies using CPAN

Depending on the version of installed Perl environment you may need to install additional modules (listed in the **Prerequisites**) from the Comprehensive Perl Archive Network (CPAN) – the repository of Perl modules. To install Perl modules use the CPAN client which should be installed by default with the Perl environment. To use the CPAN client, open the command-line interface, type the following command and press **Enter**:

```
cpan install <module name>
```

### Launching software

To launch the software simply run the following command:

```
perl vcf2caps.pl
```

### Run analysis on test data

To check if the software runs properly, perform analysis on the test data. They are available in the 'test' directory and contains the following files:
- test.vcf – contains the fragment of VCF file with identified SNPs/InDels from chromosome 1 of sugar beet;
- enzymes – contains enzyme database in GCG format.
The reference sequence of sugar beet chromosome 1 you can get from the NCBI database: 
[https://www.ncbi.nlm.nih.gov/nuccore/CM002321.2](https://www.ncbi.nlm.nih.gov/nuccore/CM002321.2)

After downloading the reference file, please confirm that the sequence name contains 'CM002321.2'.

### Getting help

Any issues/problems/comments can be posted on Github issues or you can contact me directly through the following email address: w.wesolowski at protonmail dot com.

I will try to solve all the problems you encountered using VCF2CAPS software.

### License

This software is licensed under the GNU General Public License v3.0. See LICENSE.txt for more information.
