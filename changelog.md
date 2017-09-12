# 2017-09-12 #
* Workaround for bug in MEME's XML output, introduced in MEME v4.11.x
* Add support for FIMO v4.11.4 and newer

# 2017-08-14 #
* Add citing information to readme

# 2017-08-03 #
* [web] move source code, binaries and changelog to github

# 2017-02-17 #
* [web] attach promoter positions file to results email
* [web] remove absolute paths from log file
* [web] properly remove old files from results dir

# 2016-12-23 #
* small adjustment to copyright text
* add check if (adjusted) promoter set still has at least 4 sequences if (original) set exceeded upstream cluster border

# 2016-10-26 #
* [web] add a more clear phrase to the result email, if no cluster was predicted and no error occurred
* [web] always attach log file to result email

# 2016-10-20 #
* fix problem with different input line seperators in files created on Linux vs. files created on Windows machines
* fix small typo in upstream cluster extension code (if gap length == 5)
* [web] improve handling of filenames without extensions
* [web] job submission summary: add error message if anchor gene does not occur in annotation file

# 2016-07-26 #
* fix bug in promoter sequence calculation
* improve error messages
* fix bug which causes false reports of faulty promoters in cluster predictions
* improve some regular expressions to avoid problems with meta-characters in gene names
* drop support for AspGD, JGI, GFF, … annotation file formats (caused quite some trouble)
* introduce [self-defined annotation file format](https://sbi.hki-jena.de/cassis/Help.php#Input)

# 2016-07-20 #
* [web] improve handling of unusual file names

# 2016-04-21 #
* CASSIS is now able to work with [MEME suite versions newer than v4.9.1](http://meme-suite.org/doc/download.html) (preserving backwards compatibility with older versions, like v4.9.1)
* improve handling of tainted $PATH variables

# 2016-04-01 #
* provide [SiTaR](http://www.ncbi.nlm.nih.gov/pubmed/21893518) as an alternative to MEME/FIMO, if TFBSs are already known and don't need to be predicted de novo
* better error handling
* small improvements
