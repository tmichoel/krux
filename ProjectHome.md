kruX is a toolbox for performing millions of [non-parametric ANOVA (Kruskal-Wallis)](http://en.wikipedia.org/wiki/Kruskal–Wallis_one-way_analysis_of_variance) tests at once using matrix-multiplication methods, implemented in Matlab, Python and R.

Imagine you have observations for a large number of variables which can be subdivided into a large number of possible groupings. You want to test for each variable whether or not the samples originate from the same distribution when divided according to each of the possible groupings. For instance, in [eQTL](http://en.wikipedia.org/wiki/Expression_quantitative_trait_loci) studies the number of variables can run into the ten-thousands and the number of possible groupings can be upto a few million. Although performing a single Kruskal-Wallis test in any of the popular statistics software packages may not take longer than a few thousands of a second, testing this many combinations easily takes you into months of single-CPU run times. In contrast, kruX will do the same task in not more than **40 minutes** (typical human dataset) of single-CPU run time.

Credit where credit is due: development of kruX was inspired by [Matrix-eQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) which performs the same task using [parametric ANOVA](http://en.wikipedia.org/wiki/ANOVA) tests.

The following links will get you started with kruX:

  * [Download and installation instructions](https://code.google.com/p/krux/wiki/InstallationInstructions)
  * [Matlab tutorial](https://krux.googlecode.com/svn/trunk/html/kruX_tutorial_matlab.html)

> For more information, contact <a href='http://lab.michoel.info'>Tom Michoel</a>.

**News**

  * 14 Jan 2013 - kruX paper final version published in [BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/15/11/) and also posted on the [arXiv](http://arxiv.org/abs/1307.3519).
  * 24 Oct 2013 - kruX paper revision posted on the [arXiv](http://arxiv.org/abs/1307.3519).
  * 27 Sep 2013 - major software update (Matlab version), 4x faster than previous version.
  * 12 Jul 2013 - kruX paper posted on the [arXiv](http://arxiv.org/abs/1307.3519) with more details about the algorithm and an application on a human dataset.
  * 27 May 2013 - kruX code publically available.