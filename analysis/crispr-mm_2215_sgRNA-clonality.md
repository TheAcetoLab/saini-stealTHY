---
title: "Clonal dynamics using mock CRISPR screen"
subtitle: "Mouse libray (n = 2215)"
author: "Francesc Castro-Giner"
date: 'October 31, 2025'
output: workflowr::wflow_html
editor_options:
  chunk_output_type: console
params:
  data_dir: ./data/crispr/mm_2215_sgRNA
---

<p>
<button type="button" class="btn btn-default btn-workflowr btn-workflowr-report"
  data-toggle="collapse" data-target="#workflowr-report">
  <span class="glyphicon glyphicon-list" aria-hidden="true"></span>
  workflowr
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
</button>
</p>

<div id="workflowr-report" class="collapse">
<ul class="nav nav-tabs">
  <li class="active"><a data-toggle="tab" href="#summary">Summary</a></li>
  <li><a data-toggle="tab" href="#checks">
  Checks <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  </a></li>
  <li><a data-toggle="tab" href="#versions">Past versions</a></li>
</ul>

<div class="tab-content">
<div id="summary" class="tab-pane fade in active">
  <p><strong>Last updated:</strong> 2025-10-31</p>
  <p><strong>Checks:</strong>
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  7
  <span class="glyphicon glyphicon-exclamation-sign text-danger" aria-hidden="true"></span>
  0
  </p>
  <p><strong>Knit directory:</strong>
  <code>saini-stealTHY/</code>
  <span class="glyphicon glyphicon-question-sign" aria-hidden="true"
  title="This is the local directory in which the code in this file was executed.">
  </span>
  </p>
  <p>
  This reproducible <a href="https://rmarkdown.rstudio.com">R Markdown</a>
  analysis was created with <a
  href="https://github.com/workflowr/workflowr">workflowr</a> (version
  1.7.1). The <em>Checks</em> tab describes the
  reproducibility checks that were applied when the results were created.
  The <em>Past versions</em> tab lists the development history.
  </p>
<hr>
</div>
<div id="checks" class="tab-pane fade">
  <div class="panel-group" id="workflowr-checks">
  <div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongRMarkdownfilestronguptodate">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>R Markdown file:</strong> up-to-date
</a>
</p>
</div>
<div id="strongRMarkdownfilestronguptodate" class="panel-collapse collapse">
<div class="panel-body">
  
Great! Since the R Markdown file has been committed to the Git repository, you
know the exact version of the code that produced these results.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongEnvironmentstrongempty">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Environment:</strong> empty
</a>
</p>
</div>
<div id="strongEnvironmentstrongempty" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! The global environment was empty. Objects defined in the global
environment can affect the analysis in your R Markdown file in unknown ways.
For reproduciblity it's best to always run the code in an empty environment.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongSeedstrongcodesetseed20240517code">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Seed:</strong> <code>set.seed(20240517)</code>
</a>
</p>
</div>
<div id="strongSeedstrongcodesetseed20240517code" class="panel-collapse collapse">
<div class="panel-body">
  
The command <code>set.seed(20240517)</code> was run prior to running the code in the R Markdown file.
Setting a seed ensures that any results that rely on randomness, e.g.
subsampling or permutations, are reproducible.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongSessioninformationstrongrecorded">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Session information:</strong> recorded
</a>
</p>
</div>
<div id="strongSessioninformationstrongrecorded" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! Recording the operating system, R version, and package versions is
critical for reproducibility.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongCachestrongnone">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Cache:</strong> none
</a>
</p>
</div>
<div id="strongCachestrongnone" class="panel-collapse collapse">
<div class="panel-body">
  
Nice! There were no cached chunks for this analysis, so you can be confident
that you successfully produced the results during this run.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongFilepathsstrongrelative">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>File paths:</strong> relative
</a>
</p>
</div>
<div id="strongFilepathsstrongrelative" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! Using relative paths to the files within your workflowr project
makes it easier to run your code on other machines.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongRepositoryversionstrongahrefhttpsgithubcomTheAcetoLabsainistealTHYtree455fa9121fcca51363cf639f7ee170fb108e34c1targetblank455fa91a">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Repository version:</strong> <a href="https://github.com/TheAcetoLab/saini-stealTHY/tree/455fa9121fcca51363cf639f7ee170fb108e34c1" target="_blank">455fa91</a>
</a>
</p>
</div>
<div id="strongRepositoryversionstrongahrefhttpsgithubcomTheAcetoLabsainistealTHYtree455fa9121fcca51363cf639f7ee170fb108e34c1targetblank455fa91a" class="panel-collapse collapse">
<div class="panel-body">
  
<p>
Great! You are using Git for version control. Tracking code development and
connecting the code version to the results is critical for reproducibility.
</p>

<p>
The results in this page were generated with repository version <a href="https://github.com/TheAcetoLab/saini-stealTHY/tree/455fa9121fcca51363cf639f7ee170fb108e34c1" target="_blank">455fa91</a>.
See the <em>Past versions</em> tab to see a history of the changes made to the
R Markdown and HTML files.
</p>

<p>
Note that you need to be careful to ensure that all relevant files for the
analysis have been committed to Git prior to generating the results (you can
use <code>wflow_publish</code> or <code>wflow_git_commit</code>). workflowr only
checks the R Markdown file, but you know if there are other scripts or data
files that it depends on. Below is the status of the Git repository when the
results were generated:
</p>

<pre><code>
Ignored files:
	Ignored:    .DS_Store
	Ignored:    .Rhistory
	Ignored:    .Rproj.user/
	Ignored:    analysis/.DS_Store
	Ignored:    analysis/templates/.DS_Store
	Ignored:    code/.DS_Store
	Ignored:    code/raw_data_processing/rnaseq/p27851_o32062/pipelines/
	Ignored:    configuration/.DS_Store
	Ignored:    data/.DS_Store
	Ignored:    data/crispr/
	Ignored:    data/resources/
	Ignored:    data/rnaseq/
	Ignored:    output/.DS_Store
	Ignored:    output/clinical/
	Ignored:    output/crispr/
	Ignored:    output/rnaseq/

Untracked files:
	Untracked:  analysis/index.md
	Untracked:  analysis/style.css
	Untracked:  code/R-functions/subchunkify.R

Unstaged changes:
	Modified:   .gitignore
	Modified:   analysis/_site.yml
	Deleted:    analysis/crispr-hsapiens_2180_sgRNA_r1.Rmd
	Deleted:    analysis/crispr-mm_2215_sgRNA-clonality_r1.Rmd
	Deleted:    analysis/crispr-mm_2215_sgRNA-r2oC.Rmd
	Deleted:    analysis/crispr-mm_2215_sgRNA-r2oD.Rmd
	Deleted:    analysis/crispr-mm_2215_sgRNA-r2oD_r1.Rmd
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_ko_culture--over--amhr2_control_culture_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_ko_culture--over--amhr2_control_culture_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_ko_notgfbeta--over--amhr2_ko_culture_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_ko_notgfbeta--over--amhr2_ko_culture_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_control_culture_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_control_culture_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_ko_culture_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_ko_culture_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_ko_notgfbeta_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-amhr2_overexpression--over--amhr2_ko_notgfbeta_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/figure/gse-dotplot-invivo_amhr2_ko--over--invivo_control_sgrna_GSEA_msigdb.h-1.pdf
	Deleted:    analysis/templates/docs/figure/gse-dotplot-invivo_amhr2_ko--over--invivo_control_sgrna_GSEA_msigdb.h-1.png
	Deleted:    analysis/templates/docs/file/tmp/dge.xlsx
	Deleted:    analysis/templates/docs/file/tmp/gse_gsea_msigdb.h.xlsx
	Deleted:    analysis/templates/docs/file/tmp/gse_gsego_bp.xlsx
	Deleted:    analysis/templates/docs/file/tmp/gse_gsego_mf.xlsx
	Modified:   analysis/templates/rnaseq-deg-edger.Rmd
	Modified:   code/R-functions/gse_report.r
	Modified:   configuration/rmarkdown/color_palettes.R
	Modified:   configuration/rmarkdown/ggplot_theme.R
	Modified:   update_workflowr.R

</code></pre>

<p>
Note that any generated files, e.g. HTML, png, CSS, etc., are not included in
this status report because it is ok for generated content to have uncommitted
changes.
</p>

</div>
</div>
</div>
</div>
<hr>
</div>
<div id="versions" class="tab-pane fade">
  
<p>
These are the previous versions of the repository in which changes were made
to the R Markdown (<code>analysis/crispr-mm_2215_sgRNA-clonality.Rmd</code>) and HTML (<code>docs/crispr-mm_2215_sgRNA-clonality.html</code>)
files. If you've configured a remote Git repository (see
<code>?wflow_git_remote</code>), click on the hyperlinks in the table below to
view the files as they were in that past version.
</p>
<div class="table-responsive">
<table class="table table-condensed table-hover">
<thead>
<tr>
<th>File</th>
<th>Version</th>
<th>Author</th>
<th>Date</th>
<th>Message</th>
</tr>
</thead>
<tbody>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/455fa9121fcca51363cf639f7ee170fb108e34c1/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">455fa91</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Prepare files for publication</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/crispr-mm_2215_sgRNA-clonality.html" target="_blank">a6972e9</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/b6bc77088b14bf2712a35dc41be83f6032057e27/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">b6bc770</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>update crispr analysis</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/f982e709689944cbcd3283f0040d2811c625c7d4/docs/crispr-mm_2215_sgRNA-clonality.html" target="_blank">f982e70</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/f2497227bea4d0c618e871f06ff31d6181a6e3d2/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">f249722</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>mod clonality figuress</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/0e8d1db849f537c8abb88e943b7f6fd049d838fc/docs/crispr-mm_2215_sgRNA-clonality.html" target="_blank">0e8d1db</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/748478ddf826aace4d4ea19c1a8a6bf5c7f1ba6e/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">748478d</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>mod clonality figuress</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/e4cd8b5b49c70c2dac6df9053e5582afcfcd89c5/docs/crispr-mm_2215_sgRNA-clonality.html" target="_blank">e4cd8b5</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/ab1136cc72c79ad56e6d477f4078b442e134505d/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">ab1136c</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-13</td>
<td>update clonality</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/fa3898f6c69fe480547f4cc28cc5a78fce514b67/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">fa3898f</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-12</td>
<td>updated revisions</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e513b9973e91772c074941bf5bcfb1a1eaf0380e/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">e513b99</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-08</td>
<td>update clonality figures</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/2dda789fcb78485deaf0274d2ef9be581c62e90c/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">2dda789</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-07</td>
<td>update clonality figures</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/382cde0b92a2a809b4848d543a398e0d89ffc454/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">382cde0</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-07</td>
<td>update clonality figures</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/c2a5c1c6d68eaf53b7e3abc93f74bd519068b3c2/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">c2a5c1c</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-06</td>
<td>update r1od</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/b6d9355f3511e9503b4f6da1b3b2571277a4a223/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">b6d9355</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-01</td>
<td>mod r2oD</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/5fab19d765a869233e820f2aed94595a89d5fe9b/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">5fab19d</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-28</td>
<td>clean code</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/1920712c9ede8023880e174fc87f3cd238c3d612/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">1920712</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-24</td>
<td>added F5</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/628f8ec8600fd6671cb7c7f57c8e717199ff7f75/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">628f8ec</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-24</td>
<td>added F4</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/abf695c6a645a3aa42ea9e2fb622b5f22e723b11/analysis/crispr-mm_2215_sgRNA-clonality.Rmd" target="_blank">abf695c</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-23</td>
<td>done f2 anf f3</td>
</tr>
</tbody>
</table>
</div>

<hr>
</div>
</div>
</div>






## Load libraries, additional functions and data

Setup environment

``` r
knitr::opts_chunk$set(results='asis', echo=TRUE, message=FALSE, warning=FALSE, error=FALSE, fig.align = 'center', fig.width = 3.5, fig.asp = 0.618, dpi = 600, dev = c("png", "pdf"), fig.showtext = FALSE, engine.opts = list(bash = "-l"))

options(stringsAsFactors = FALSE)

use_seed <- 1100101

set.seed(use_seed)
```

Load packages

``` r
library(tidyverse)
library(knitr)
library(foreach)
library(magrittr)
library(DT)
library(kableExtra)
library(diptest)
library(SummarizedExperiment)

library(ggridges)
library(ggh4x)
library(patchwork)
library(colorblindr)
library(ggbeeswarm)
library(ggpubr)
```

Load ggplot theme

``` r
source("./configuration/rmarkdown/ggplot_theme.R")
```

Load ggplot theme

``` r
source("./configuration/rmarkdown/color_palettes.R")
```

Load Summarized Experiment object

``` r
se <- readRDS(file.path(params$data_dir, 'se.rds'))
```

Filter samples used for this analysis

``` r
use_cols <- colData(se) %>% 
  data.frame %>% 
  filter(effect_type_short == 'noKO') %>% 
  filter(vector_short != "Thy1_EF1alpha_dCas9") %>% 
  rownames

se <- se[,use_cols]
```

Data wrangling

``` r
chunk_colData <- colData(se) %>% data.frame %>%
  mutate(
    vector_short = ifelse(vector_short == 'Puro_EF1alpha_SpCas9', 
                          'Cas9+Puro', vector_short),
    vector_short = ifelse(vector_short == 'Puro+ dCas9', 
                          'Cas9+Puro', vector_short),
    vector_short = factor(vector_short, levels = c('Thy1', 'Puro', 'Cas9+Puro')),
    mouse_model_short = case_match(mouse_model_global,
                                   'immunodeficient' ~ 'NSG',
                                   'immunocompetent' ~ 'Syngeneic',
                                   'C57BL6' ~ 'Syngeneic'
                                   ),
    mouse_model_short = factor(mouse_model_short, levels = c('NSG', 'Syngeneic')),
    vector_mhost = paste0(vector_short, ', ', mouse_model_short),
    model_vector_mhost = paste0(donor,', ', vector_short, ', ', mouse_model_short),
    sample_type = ifelse(cancer_type %in% c('Breast', 'CRC'),
                         case_match(sample_type, 
                             'whole_tumor' ~ 'Primary tumor', 
                             'whole_lung' ~ 'Lung mets',
                             'bulk' ~ 'Bulk',
                             'cell_culture' ~ 'Cell culture'),
                         sample_type
                         ),
        sample_type = ifelse(cancer_type %in% c('NSCLC'),
                         case_match(sample_type, 
                             'met_liver' ~ 'Liver mets', 
                             'whole_lung' ~ 'Primary tumor',
                             'bulk' ~ 'Bulk',
                             'cell_culture' ~ 'Cell culture'),
                         sample_type
                         ),
    
    sample_type = factor(sample_type, levels = c('Primary tumor', 'Lung mets', 'Liver mets', 'Bulk', 'Cell culture'))
    
  )

colData(se) <- chunk_colData %>% DataFrame
```

## Configure analyses

Define comparisons for modality and Kolmogorov-Smirnov test

``` r
x <- colData(se) %>% data.frame

group_list_modality <- list(
  `Syngeneic Thy1` = x %>% 
    filter(vector_mmodel == '4T1-noKO-Thy1-BALB' & sample_type == 'Primary tumor') %>% 
    pull(sample_alias),
  
  `Syngeneic Puro` = x %>% 
    filter(vector_mmodel %in% c('4T1-noKO-Puro-BALB', 'MVT1-noKO-Puro-FVB', 'Py2T-noKO-Puro-FVB') & sample_type == 'Primary tumor') %>%
    pull(sample_alias),
  
  `Syngeneic dCas9 + Puro` = x %>% 
    filter(vector_mmodel == '4T1-noKO-Puro_EF1alpha_dCas9-BALB' & sample_type == 'Primary tumor') %>% 
    pull(sample_alias),
  
  `NSG Thy1` =  x %>% 
    filter(vector_mmodel %in% c('4T1-noKO-Thy1-NSG', 'MVT1-noKO-Thy1-NSG', 'Py2T-noKO-Thy1-NSG') & sample_type == 'Primary tumor') %>%
    pull(sample_alias),
  
  `NSG Puro` = x %>% 
    filter(vector_mmodel %in% c('4T1-noKO-Puro-NSG', 'MVT1-noKO-Puro-NSG', 'Py2T-noKO-Puro-NSG') & sample_type == 'Primary tumor') %>%
    pull(sample_alias),
  
  `NSG dCas9 + Puro` = x %>% 
    filter(vector_mmodel == '4T1-noKO-Puro_EF1alpha_dCas9-NSG' & sample_type == 'Primary tumor') %>% 
    pull(sample_alias)
)

comp_list_modality <- list(
  `Syngeneic Thy1 vs. dCas9 + Puro` = list(
    `dCas9 + Puro` = group_list_modality$`Syngeneic dCas9 + Puro`,
    `Thy1` = group_list_modality$`Syngeneic Thy1`
  ),
  `Syngeneic Thy1 vs. Puro` = list(
    `Puro` = group_list_modality$`Syngeneic Puro`,
    `Thy1` = group_list_modality$`Syngeneic Thy1`
  )
)

group_df_modality <- foreach(i = names(group_list_modality), .combine = rbind) %do% {
  data.frame(
    sample_alias = group_list_modality[[i]],
    modality_group = i
  )
}
```


## sgRNA distribution in Primary tumor

``` r
use_sample_type <- 'Primary tumor'
use_colData <- colData(se) %>% data.frame %>%
  filter(sample_type == use_sample_type)
```

### Figure 2B: sgRNA distribution in mice transplanted with 4T1 cells

Density plots showing unique sgRNAs in 4T1 tumors (n=5); frequency = read counts per million (cpm+1).


``` r
selected_panels <- list(
  `Thy1, NSG` = c('4T_TumoT_N4', '4T_TumoT_N5', '4T_TumoT_N2'),
  `Thy1, Syngeneic` = c('4T_TumoT_B1', '4T_TumoT_B3', '4T_TumoT_B4'),
 
  `Puro, NSG` = c('4T_TumoP_N5', '4T_TumoP_N4', '4T_TumoP_N2'),
  `Puro, Syngeneic` = c('4T_TumoP_B3', '4T_TumoP_B2', '4T_TumoP_B1'),
  
  `dCas9+Puro, NSG` = c('PlatePositionE1', 'PlatePositionC1', 'PlatePositionA2'),
  `dCas9+Puro, Syngeneic` = c('PlatePositionE2', 'PlatePositionG2', 'PlatePositionC2')
  
)

selected_panels_df <- foreach(i = names(selected_panels), .combine = rbind) %do% {
  data.frame(group = i, sample_alias = selected_panels[[i]])
}

use_samples <- unlist(selected_panels)
use_assay <- assay(se[,use_samples], 'cpm') %>% 
  as.data.frame(check.names = F) %>% 
  rownames_to_column('guide') %>% 
  pivot_longer(-guide, names_to = 'sample_alias', values_to = 'cpm') %>% 
  left_join(use_colData) %>% 
  left_join(selected_panels_df) %>% 
  mutate(
    sample_alias = factor(sample_alias, levels = use_samples),
    group = factor(group, levels = names(selected_panels))
  )

use_assay %>% 
  ggplot(aes(x=cpm + 1, y = sample_alias, 
             fill = vector_short, color = vector_short, 
             height = after_stat(density))) +
  geom_density_ridges(size = one_pt/4, 
                      alpha = 0.8, 
                      scale = 4,
                      rel_min_height = 0.001 # set the `rel_min_height` argument to remove tails,
                      ) +
  scale_fill_manual(values = palette_vector) +
  scale_color_manual(values = palette_vector_line) +
  scale_x_log10(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = c(0.01, 0)) +
  labs(y = '') +
  guides(fill = 'none', color = 'none') +
  facet_wrap2(vars(group), ncol = 2,
                     scales = 'free_y',
                     axes = "x") +
  theme_ridges(font_size = 3, grid = F) +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_line(linewidth = one_pt/4, color = 'black'),
    axis.ticks.x = element_line(linewidth = one_pt/4, color = 'black'),
    axis.text.y = element_blank(),
    strip.text.x = element_text(size=4, hjust = 0)
  )
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/test-dist-tumor-distribution-ridges_plot-4t1-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-test-dist-tumor-distribution-ridges_plot-4t1-1">
  Past versions of test-dist-tumor-distribution-ridges_plot-4t1-1.png
  </button>
  </p>

  <div id="fig-test-dist-tumor-distribution-ridges_plot-4t1-1" class="collapse">
  <div class="table-responsive">
  <table class="table table-condensed table-hover">
  <thead>
  <tr>
  <th>Version</th>
  <th>Author</th>
  <th>Date</th>
  </tr>
  </thead>
  <tbody>
  <tr>
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e4cd8b5b49c70c2dac6df9053e5582afcfcd89c5/docs/figure/crispr-mm_2215_sgRNA-clonality.Rmd/test-dist-tumor-distribution-ridges_plot-4t1-1.png" target="_blank">e4cd8b5</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-05-13</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


### Supplementary Figure 2B: sgRNA distribution in mice transplanted with Py2T and MVT1 cells

Density plots showing the distribution of unique sgRNAs retrieved in primary tumors from the indicated mice transplanted with Py2T and MVT1 cells, previously transduced as indicated; frequency=read counts per million (cpm+1); each plot is representative of a group of n=5. 


``` r
selected_panels <- list(
  `Py2T, Thy1, NSG` = c('PY_TumoT_N5', 'PY_TumoT_N4', 'PY_TumoT_N3'),
  `Py2T, Thy1, Syngeneic` = c('PY_TumoT_F3', 'PY_TumoT_F4', 'PY_TumoT_F5'), 
  
  `Py2T, Puro, NSG` = c('PY_TumoP_N1', 'PY_TumoP_N2', 'PY_TumoP_N3'),
  `Py2T, Puro, Syngeneic` = c('PY_TumoP_F1', 'PY_TumoP_F5', 'PY_TumoP_F2'),
  
  `MVT1, Thy1, NSG` = c('MVT_TumoT_N3', 'MVT_TumoT_N4', 'MVT_TumoT_N5'),
  `MVT1, Thy1, Syngeneic` = c('MVT_TumoT_F1', 'MVT_TumoT_F4', 'MVT_TumoT_F3'),
  
  `MVT1, Puro, NSG` = c('MVT_TumoP_N3', 'MVT_TumoP_N4', 'MVT_TumoP_N5'),
  `MVT1, Puro, Syngeneic` = c('MVT_TumoP_F3', 'MVT_TumoP_F4', 'MVT_TumoP_F5')
)

selected_panels_df <- foreach(i = names(selected_panels), .combine = rbind) %do% {
  data.frame(group = i, sample_alias = selected_panels[[i]])
}

use_samples <- unlist(selected_panels)
use_assay <- assay(se[,use_samples], 'cpm') %>% 
  as.data.frame(check.names = F) %>% 
  rownames_to_column('guide') %>% 
  pivot_longer(-guide, names_to = 'sample_alias', values_to = 'cpm') %>% 
  left_join(use_colData) %>% 
  left_join(selected_panels_df) %>% 
  mutate(
    sample_alias = factor(sample_alias, levels = use_samples),
    group = factor(group, levels = names(selected_panels))
  )

use_assay %>% 
  ggplot(aes(x=cpm + 1, y = sample_alias, 
             fill = vector_short, 
             color = vector_short, 
             height = after_stat(density))) +
  geom_density_ridges(size = one_pt/4, 
                      alpha = 0.8, 
                      scale = 4,
                      rel_min_height = 0.001 # set the `rel_min_height` argument to remove tails,
                      ) +
  scale_fill_manual(values = palette_vector) +
  scale_color_manual(values = palette_vector_line) +
  scale_x_log10(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = c(0.01, 0)) +
  labs(y = '') +
  guides(fill = 'none', color = 'none') +
  facet_wrap2(vars(group), ncol = 2,
                     scales = 'free_y',
                     axes = "x") +
  theme_ridges(font_size = 3, grid = F) +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_line(linewidth = one_pt/4, color = 'black'),
    axis.ticks.x = element_line(linewidth = one_pt/4, color = 'black'),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill =NULL,
                                color = "black",
                                linewidth = 2*one_pt)
  )
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/test-dist-tumor-distribution-ridges_plot-py2t-mvt1-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-test-dist-tumor-distribution-ridges_plot-py2t-mvt1-1">
  Past versions of test-dist-tumor-distribution-ridges_plot-py2t-mvt1-1.png
  </button>
  </p>

  <div id="fig-test-dist-tumor-distribution-ridges_plot-py2t-mvt1-1" class="collapse">
  <div class="table-responsive">
  <table class="table table-condensed table-hover">
  <thead>
  <tr>
  <th>Version</th>
  <th>Author</th>
  <th>Date</th>
  </tr>
  </thead>
  <tbody>
  <tr>
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e4cd8b5b49c70c2dac6df9053e5582afcfcd89c5/docs/figure/crispr-mm_2215_sgRNA-clonality.Rmd/test-dist-tumor-distribution-ridges_plot-py2t-mvt1-1.png" target="_blank">e4cd8b5</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-05-13</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


## Kolmogorov-Smirnov test and ECDF plots 

Run Kolmogorov-Smirnov test to compare distributions merged by group and create Empirical Cumulative Density Function (ECDF) by comparison.

``` r
i <- names(comp_list_modality)[1]
ks_test_df <- data.frame()
ecdf_plots <- list()
for(i in names(comp_list_modality)) {
  use_samples <- comp_list_modality[[i]]
  use_assay_c <- assay(se[,use_samples[[1]]], 'cpm') %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column('guide') %>% 
    pivot_longer(-guide, values_to = 'cpm',names_to = 'sample_alias') %>% 
    mutate(group = 'case', condition = names(use_samples)[1])
  use_assay_r <- assay(se[,use_samples[[2]]], 'cpm') %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column('guide') %>% 
    pivot_longer(-guide, values_to = 'cpm',names_to = 'sample_alias') %>% 
    mutate(group = 'referene', condition = names(use_samples)[2])
  use_assay <- rbind(use_assay_c, use_assay_r)
  
  # Perform the Kolmogorov-Smirnov test
  ks.res <- ks.test(use_assay_c$cpm, use_assay_r$cpm)
  ks_test_df <- rbind(
    ks_test_df,
      data.frame(
        comparison = i,
        case = names(use_samples)[1],
        reference = names(use_samples)[2],
        D = ks.res$statistic,
        P = ks.res$p.value
      )
    )

  # Create ECDFs plot Empirical Cumulative Density Function
  # use_colors <- c('red', 'blue') %>% set_names(names(use_samples))
  ecdf_plots[[i]] <- use_assay %>% 
    ggplot(aes(x = log(1+cpm), color = condition, group = sample_alias)) +
    stat_ecdf(geom = 'step', alpha = 1, linewidth = one_pt/4) +
    scale_color_manual(values = palette_vector) +
    labs(x = "log(1+cpm)",
         color = '',
         sub = i
         ) +
    theme(
      text = element_text(size=2),
      axis.text.y =  element_text(size=3),
      axis.text.x =  element_text(size=3),
      legend.text =  element_text(size=2),
      legend.title =  element_text(size=2),
      plot.subtitle = element_text(size=3),
      plot.caption =  element_text(size=2)
      
    )

}
```

### Supplementary Figure 2C left : K-S test table 

Table depicting the results of the Kolmogorov-Smirnov test for comparing sgRNA distributions from the datasets indicated, encompassing three mammary carcinoma models (MVT1, 4T1, Py2T). The datasets are obtained from primary tumors shown in Figure 2B and Supplementary Figure 2B.


``` r
ks_test_df %>% 
  dplyr::select(-case, -reference) %>% 
  mutate(P = format.pval(P)) %>% 
  datatable(., 
            rownames = FALSE, 
            filter = 'top', 
            caption = 'Results of Kolmogorov-Smirnov test',
            extensions = 'Buttons', 
            options = list(
              dom = 'Blfrtip',
              buttons = c('csv', 'excel'),
              columnDefs = list(
                list(width = 120, targets = 1)
              )
            )
            ) %>% 
    formatRound(columns = c('D'), digits = 3)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-0e41329c1963262c826a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0e41329c1963262c826a">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.375289691497238\" data-max=\"0.732731376975173\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>Results of Kolmogorov-Smirnov test<\/caption>","data":[["Syngeneic Thy1 vs. dCas9 + Puro","Syngeneic Thy1 vs. Puro"],[0.7327313769751728,0.3752896914972385],["&lt; 2.22e-16","&lt; 2.22e-16"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>comparison<\/th>\n      <th>D<\/th>\n      <th>P<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"width":120,"targets":1},{"className":"dt-right","targets":1},{"name":"comparison","targets":0},{"name":"D","targets":1},{"name":"P","targets":2}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":["options.columnDefs.0.render"],"jsHooks":[]}</script>
```

### Supplementary Figure 2C right : Empirical Cumulative Density Function (ECDF) plots

Cumulative curves depicting the results of the Kolmogorov-Smirnov test for comparing sgRNA distributions from the datasets indicated, encompassing three mammary carcinoma models (MVT1, 4T1, Py2T). The datasets are obtained from primary tumors shown in Figure 2B and Supplementary Figure 2B.

Empirical Cumulative Density Function (ECDF) plot by comparison

``` r
use_index <- c("Syngeneic Thy1 vs. Puro", "Syngeneic Thy1 vs. dCas9 + Puro")
plot_grid(plotlist = ecdf_plots[use_index], ncol = 2)
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/modality-comparisons-ks-ecdf-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-modality-comparisons-ks-ecdf-1">
  Past versions of modality-comparisons-ks-ecdf-1.png
  </button>
  </p>

  <div id="fig-modality-comparisons-ks-ecdf-1" class="collapse">
  <div class="table-responsive">
  <table class="table table-condensed table-hover">
  <thead>
  <tr>
  <th>Version</th>
  <th>Author</th>
  <th>Date</th>
  </tr>
  </thead>
  <tbody>
  <tr>
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e4cd8b5b49c70c2dac6df9053e5582afcfcd89c5/docs/figure/crispr-mm_2215_sgRNA-clonality.Rmd/modality-comparisons-ks-ecdf-1.png" target="_blank">e4cd8b5</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-05-13</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

## Figure 2C: Representative Müller plots of 4T1 in syngeneic recipients

Define selected comparisons for muller plots

``` r
comp_list_muller_syngeneic_representative <- list(
 `Thy1` = c('4T_TumoT_B3', '4T_LungsT_B3'),
 `Puro` = c('4T_TumoP_B1', '4T_LungsP_B1'),
 `dCas9 + Puro` = c('PlatePositionA3', 'PlatePositionB3')
)
```

Generate muller plots

``` r
cpm_threshold <- 10000
comp_list <- comp_list_muller_syngeneic_representative

use_df <- assay(se, 'cpm') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('guide') %>% 
  mutate(
    gene = rowData(se)$Gene,
    gene = ifelse(grepl('Non_Target', gene), 'Non_Target', gene)
  ) %>%
  pivot_longer(-c(guide, gene), names_to = 'sample_alias', values_to = 'cpm')

j <-  names(comp_list)[1]
muller_plots <- foreach(j = names(comp_list)) %do% {
  i <- comp_list[[j]]
  x <- use_df %>% 
    filter(sample_alias %in% i) %>% 
    mutate(
      sample_type = ifelse(sample_alias == i[1], 'Primary tumor', 'Lungs'),
      sample_type = factor(sample_type, levels = c('Primary tumor', 'Lungs'))
    )
  
  # Select guides with cpm > cpm_threshold
  x_guides_keep <- x %>% 
    group_by(sample_type, guide) %>% 
    summarise(cpm = max(cpm)) %>% 
    filter(cpm > cpm_threshold) %>% 
    pull(guide)
  
  # Group guides with cpm < cpm_threshold into the low abundant category, and sum counts for low abundant
  x %<>% 
    mutate(
      guide = ifelse(guide %in%  x_guides_keep, guide, 'low abundant')
    ) %>% 
    group_by(sample_type, guide) %>%
    summarise(cpm = sum(cpm)) %>% 
    ungroup() %>%
    mutate(
      guide = fct_reorder(guide, cpm),
      guide = relevel(guide, ref = 'low abundant')
    )
  
  use_cols <- c(
    'low abundant' = 'grey80',
    colorRampPalette(rev(palette_OkabeIto))(nlevels(x$guide) - 1) %>% set_names(levels(x$guide)[-1])
  )
  
  x %>% 
    ggplot( aes(x = sample_type, y = cpm/10000, group = guide, fill = guide)) + 
    geom_area(colour = alpha("white", 0.1), linewidth = 0.08, alpha = 0.8) +
    scale_fill_manual(values = use_cols, guide = guide_legend(ncol = 3)) +
    labs(
      x = '',
      y = expression("CPM x 10"^-4),
      title = j,
      fill = ''
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1000000/10000)) +
    guides(fill = FALSE) +
    theme(
      panel.background = element_rect(fill = "grey90"),
      # plot.margin = margin(0.5, 1.5, 0.5, 0, "cm"),
      plot.margin = margin(0, 0.5, 0, 0, "cm"),
      axis.text =  element_text(size=3),
      axis.title =  element_text(size=3),
      plot.title = element_text(size=3, hjust = 0.5)
    )

}
names(muller_plots) <- names(comp_list)
muller_syngeneic_representative <- muller_plots
```


Muller plots comparing sgRNAs frequency in matched primary tumors versus lung metastases from syngeneic 4T1 transplants (n=5; representative shown). 


``` r
wrap_plots(muller_plots, ncol = 1)
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-syngeneic-representative-plots-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-muller-syngeneic-representative-plots-1">
  Past versions of muller-syngeneic-representative-plots-1.png
  </button>
  </p>

  <div id="fig-muller-syngeneic-representative-plots-1" class="collapse">
  <div class="table-responsive">
  <table class="table table-condensed table-hover">
  <thead>
  <tr>
  <th>Version</th>
  <th>Author</th>
  <th>Date</th>
  </tr>
  </thead>
  <tbody>
  <tr>
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e4cd8b5b49c70c2dac6df9053e5582afcfcd89c5/docs/figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-syngeneic-representative-plots-1.png" target="_blank">e4cd8b5</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-05-13</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


## Hartigans' dip test
Hartigans' dip test for unimodality. If the p-value is less than 0.05, there's enough evidence to claim that the data are not unimodal (multimodal, at least bimodal).


``` r
use_assay <- assay(se, 'cpm')
i <- colnames(use_assay)[1]
modality_test_df <- foreach(i= colnames(use_assay), .combine = rbind) %do%{
  x <- use_assay[,i]
  dip.res <- dip.test(x, simulate.p.value = FALSE, B = 2000)
  data.frame(
    sample_alias = i,
    dip.D = dip.res$statistic,
    dip.p = dip.res$p.value,
    dip.result = ifelse(dip.res$p.value < 0.05, 'multimodal', 'unimodal')
    )
}

modality_test_df %<>% 
  left_join(group_df_modality) %>% 
  left_join(colData(se) %>% data.frame)
```

### Dip test table for Syngeneic models

``` r
modality_test_df %>% 
  filter(!is.na(modality_group)) %>% 
  filter(mouse_model_global == 'immunocompetent') %>% 
  mutate(
    modality_group = factor(modality_group, levels = names(group_list_modality))
  ) %>%
  group_by(modality_group) %>%
  summarise(
    `Median D` = median(dip.D) %>% round(3),
    n = n(),
    `Multimodal samples (n)` = sum(dip.result == 'multimodal'),
    `Unimodal samples (n)` = sum(dip.result == 'unimodal'),
    `Multimodal samples (%)` = round(100*`Multimodal samples (n)`/n, 1),
    `Unimodal samples (%)` = round(100*`Unimodal samples (n)`/n, 1)
  ) %>% 
  data.frame(check.names = F) %>% 
  kbl(caption = 'Dip test results for Syngeneic models') %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>Dip test results for Syngeneic models</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> modality_group </th>
   <th style="text-align:right;"> Median D </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> Multimodal samples (n) </th>
   <th style="text-align:right;"> Unimodal samples (n) </th>
   <th style="text-align:right;"> Multimodal samples (%) </th>
   <th style="text-align:right;"> Unimodal samples (%) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Syngeneic Thy1 </td>
   <td style="text-align:right;"> 0.005 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 20.0 </td>
   <td style="text-align:right;"> 80.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Syngeneic Puro </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 93.3 </td>
   <td style="text-align:right;"> 6.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Syngeneic dCas9 + Puro </td>
   <td style="text-align:right;"> 0.077 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 100.0 </td>
   <td style="text-align:right;"> 0.0 </td>
  </tr>
</tbody>
</table>

### Dip test table for NSG models

``` r
modality_test_df %>% 
  filter(!is.na(modality_group)) %>% 
  filter(mouse_model_global == 'immunodeficient') %>% 
  mutate(
    modality_group = factor(modality_group, levels = names(group_list_modality))
  ) %>%
  group_by(modality_group) %>%
  summarise(
    `Median D` = median(dip.D) %>% round(3),
    n = n(),
    `Multimodal samples (n)` = sum(dip.result == 'multimodal'),
    `Unimodal samples (n)` = sum(dip.result == 'unimodal'),
    `Multimodal samples (%)` = round(100*`Multimodal samples (n)`/n, 1),
    `Unimodal samples (%)` = round(100*`Unimodal samples (n)`/n, 1)
  ) %>% 
  data.frame(check.names = F) %>% 
  kbl(caption = 'Dip test results for NSG models') %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>Dip test results for NSG models</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> modality_group </th>
   <th style="text-align:right;"> Median D </th>
   <th style="text-align:right;"> n </th>
   <th style="text-align:right;"> Multimodal samples (n) </th>
   <th style="text-align:right;"> Unimodal samples (n) </th>
   <th style="text-align:right;"> Multimodal samples (%) </th>
   <th style="text-align:right;"> Unimodal samples (%) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NSG Thy1 </td>
   <td style="text-align:right;"> 0.005 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 13.3 </td>
   <td style="text-align:right;"> 86.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSG Puro </td>
   <td style="text-align:right;"> 0.005 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0 </td>
   <td style="text-align:right;"> 100.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NSG dCas9 + Puro </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.0 </td>
   <td style="text-align:right;"> 100.0 </td>
  </tr>
</tbody>
</table>



## Supplementary Figure 4C

Muller plots depicting the relative frequency of the predominant sgRNAs in matched primary tumors (mammary fat pad) versus synchronous metastasis samples (lungs) from syngeneic recipients transplanted with the indicated cell model, which received StealTHY KO of the mouse interactome library; n=5 each (representative samples are shown).


``` r
se <- readRDS(file.path(params$data_dir, 'se.rds'))
comp_list <- list(
 `Condition KO Thy1_all immunocompetent TUMORS--vs--Condition KO Thy1_all immunocompetent Lungs mets`=c('CRISPRKO_Clone_14', 'CRISPRKO_Clone_19')
)

cpm_threshold <- 10000

use_df <- assay(se, 'cpm') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('guide') %>% 
  mutate(
    gene = rowData(se)$Gene,
    gene = ifelse(grepl('Non_Target', gene), 'Non_Target', gene)
  ) %>%
  pivot_longer(-c(guide, gene), names_to = 'sample_alias', values_to = 'cpm')

j <-  names(comp_list)[1]
muller_plots <- foreach(j = names(comp_list)) %do% {
  i <- comp_list[[j]]
  x <- use_df %>% 
    filter(sample_alias %in% i) %>% 
    mutate(
      sample_type = ifelse(sample_alias == i[1], 'Primary tumor', 'Lungs'),
      sample_type = factor(sample_type, levels = c('Primary tumor', 'Lungs'))
    )
  
  # Select guides with cpm > cpm_threshold
  x_guides_keep <- x %>% 
    group_by(sample_type, guide) %>% 
    summarise(cpm = max(cpm)) %>% 
    filter(cpm > cpm_threshold) %>% 
    pull(guide)
  
  # Group guides with cpm < cpm_threshold into the low abundant category, and sum counts for low abundant
  x %<>% 
    mutate(
      guide = ifelse(guide %in%  x_guides_keep, guide, 'low abundant')
    ) %>% 
    group_by(sample_type, guide) %>%
    summarise(cpm = sum(cpm)) %>% 
    ungroup() %>%
    mutate(
      guide = fct_reorder(guide, cpm),
      guide = relevel(guide, ref = 'low abundant')
    )
  
  use_cols <- c(
    'low abundant' = 'grey80',
    colorRampPalette(rev(palette_OkabeIto))(nlevels(x$guide) - 1) %>% set_names(levels(x$guide)[-1])
  )
  
  x %>% 
    ggplot( aes(x = sample_type, y = cpm/10000, group = guide, fill = guide)) + 
    geom_area(colour = alpha("white", 0.1), linewidth = 0.08, alpha = 0.8) +
    scale_fill_manual(values = use_cols, guide = guide_legend(ncol = 3)) +
    labs(
      x = '',
      y = expression("CPM x 10"^-4),
      title = j,
      fill = ''
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1000000/10000)) +
    # guides(fill = FALSE) +
    theme(
        panel.background = element_rect(fill = "grey90"),
        plot.margin = margin(0, 0.5, 0, 0, "cm"),
        axis.text =  element_text(size=3),
        axis.title =  element_text(size=3),
        plot.title = element_text(size=3, hjust = 0.5),
        legend.text = element_text(size=3),
        legend.title = element_text(size=3),
        legend.position = 'bottom',
        legend.key.size = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0.1, 'mm'),
        legend.key.spacing.x = unit(0.1, 'mm')
    )

}

print(muller_plots)
```

[[1]]
<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sf4-c-1.png" style="display: block; margin: auto;" />




## Supplementary Figure 2D Lower Panels

Muller plots depicting the relative frequency of unique sgRNAs in matched primary tumors (mammary or lung orthotopic) versus synchronous metastasis samples (lung metastasis for mammary tumors, liver metastasis for lung tumors) from syngeneic and NSG recipients transplanted with the indicated models previously transduced as indicated; n=5 per group (representative samples are shown).


``` r
se <- readRDS(file.path(params$data_dir, 'se.rds'))
comp_list <- list(
  `Condition Puro_all immunoDEFICIENT TUMORS--vs--Condition Puro_all immunoDEFICIENT Lungs mets` = list(
    MVT1 = c('MVT_TumoP_N1', 'MVT_LungsP_N1'), #muller-plots-67.png
    Py2T = c('PY_TumoP_N1', 'PY_LungsP_N1')#muller-plots-77.png
    ),
  `Condition Puro_all immunocompetent TUMORS--vs--Condition Puro_all immunocompetent Lungs mets` = list(
    MVT1 = c('MVT_TumoP_F5', 'MVT_LungsP_F5'), #muller-plots-19.png
    Py2T = c('PY_TumoP_F5', 'PY_LungsP_F5') #muller-plots-25.png
    ),
  `Condition Thy1_all immunoDEFICIENT TUMORS--vs--Condition Thy1_all immunoDEFICIENT Lungs mets` = list(
    MVT1 = c('MVT_TumoT_N1', 'MVT_LungsT_N1'), #muller-plots-127.png
    Py2T = c('PY_TumoT_N1', 'PY_LungsT_N1') #muller-plots-137.png
  ),
  `Condition Thy_all immunocompetent TUMORS--vs--Condition Thy_all immunocompetent Lungs mets` = list(
    MVT1 = c('MVT_TumoT_F2', 'MVT_LungsT_F2'), #muller-plots-39.png
    Py2T = c('PY_TumoT_F3', 'PY_LungsT_F3') #muller-plots-51.png
  )
  
)

cpm_threshold <- 10000


use_df <- assay(se, 'cpm') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('guide') %>% 
  mutate(
    gene = rowData(se)$Gene,
    gene = ifelse(grepl('Non_Target', gene), 'Non_Target', gene)
  ) %>%
  pivot_longer(-c(guide, gene), names_to = 'sample_alias', values_to = 'cpm')

i <-  names(comp_list)[1]
muller_plots <- foreach(i = names(comp_list)) %do% {
  j <- comp_list[[i]][[1]]
  res <- foreach(j = comp_list[[i]]) %do% {
    x <- use_df %>% 
      filter(sample_alias %in% j) %>% 
      mutate(sample_alias = factor(sample_alias, j))
    
    # Select guides with cpm > cpm_threshold
    x_guides_keep <- x %>% 
      group_by(sample_alias, guide) %>% 
      summarise(cpm = max(cpm)) %>% 
      filter(cpm > cpm_threshold) %>% 
      pull(guide)
    
    # Group guides with cpm < cpm_threshold into the low abundant category, and sum counts for low abundant
    x %<>% 
      mutate(
         guide = ifelse(guide %in%  x_guides_keep, guide, 'low abundant')
      ) %>% 
      group_by(sample_alias, guide) %>%
      summarise(cpm = sum(cpm)) %>% 
      ungroup() %>%
      mutate(
        guide = gsub('_guide', ' sg', guide),
        guide = fct_reorder(guide, cpm),
        guide = relevel(guide, ref = 'low abundant')
      )
    
    use_cols <- c(
      'low abundant' = 'grey80',
      colorRampPalette(rev(palette_OkabeIto))(nlevels(x$guide) - 1) %>% set_names(levels(x$guide)[-1])
    )
    
    x %>% 
      ggplot( aes(x = sample_alias, y = cpm, group = guide, fill = guide)) + 
      geom_area(colour = alpha("white", 0.1), linewidth = 0.08, alpha = 0.8) +
      scale_fill_manual(values = use_cols, guide = guide_legend(ncol = 3)) +
      labs(
        x = '',
        y = 'cpm',
        fill = ''
      ) +
      theme(
        panel.background = element_rect(fill = "grey90"),
        plot.margin = margin(0.5, 1.5, 0.5, 0, "cm")
        ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,1000000))
  }
  names(res) <- names(comp_list[[i]])
  res
}
names(muller_plots) <- names(comp_list)
```


### Supplementary Figure 2D (left lower)

``` r
plot_list <- list(
  muller_plots[['Condition Puro_all immunoDEFICIENT TUMORS--vs--Condition Puro_all immunoDEFICIENT Lungs mets']][['MVT1']] + guides(fill = 'none'),
  muller_plots[['Condition Puro_all immunocompetent TUMORS--vs--Condition Puro_all immunocompetent Lungs mets']][['MVT1']]+  guides(fill = 'none'),
  muller_plots[['Condition Puro_all immunoDEFICIENT TUMORS--vs--Condition Puro_all immunoDEFICIENT Lungs mets']][['Py2T']] + guides(fill = 'none'),
  muller_plots[['Condition Puro_all immunocompetent TUMORS--vs--Condition Puro_all immunocompetent Lungs mets']][['Py2T']]+  guides(fill = 'none')
)

wrap_plots(plot_list, ncol = 2, byrow = TRUE) %>% print
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig2-d-left-lower-1.png" style="display: block; margin: auto;" />

### Supplementary Figure 2D (right lower)

``` r
plot_list <- list(
  muller_plots[['Condition Thy1_all immunoDEFICIENT TUMORS--vs--Condition Thy1_all immunoDEFICIENT Lungs mets']][['MVT1']] + guides(fill = 'none'),
  muller_plots[['Condition Thy_all immunocompetent TUMORS--vs--Condition Thy_all immunocompetent Lungs mets']][['MVT1']]+  guides(fill = 'none'),
  muller_plots[['Condition Thy1_all immunoDEFICIENT TUMORS--vs--Condition Thy1_all immunoDEFICIENT Lungs mets']][['Py2T']] + guides(fill = 'none'),
  muller_plots[['Condition Thy_all immunocompetent TUMORS--vs--Condition Thy_all immunocompetent Lungs mets']][['Py2T']]+  guides(fill = 'none')
)

wrap_plots(plot_list, ncol = 2, byrow = TRUE) %>% print
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig2-d-right-lower-1.png" style="display: block; margin: auto;" />


## Other Müller plots

Define selected comparisons for muller plots

``` r
x <- colData(se) %>% data.frame

comp_list_muller <- list(
  `Comparison 1: CMT167 Syngeneic Thy1.1` = list(
    `CMT167 Syngeneic Thy1.1 1` = c('C167_StealTHY_Tum_1',	'C167_StealTHY_Met_1'),
    `CMT167 Syngeneic Thy1.1 3` = c('C167_StealTHY_Tum_3',	'C167_StealTHY_Met_3'),
    `CMT167 Syngeneic Thy1.1 4` = c('C167_StealTHY_Tum_4',	'C167_StealTHY_Met_4'),
    `CMT167 Syngeneic Thy1.1 5` = c('C167_StealTHY_Tum_5',	'C167_StealTHY_Met_5')
  ),
  `Comparison 2: LLC1 Syngeneic Thy1.1` = list( 
   `LLC1 Syngeneic Thy1.1 1` =  c('LLC1_StealTHY_Tum_1',	'LLC1_StealTHY_Met_1'),
   `LLC1 Syngeneic Thy1.1 2` =  c('LLC1_StealTHY_Tum_2',	'LLC1_StealTHY_Met_2'),
   `LLC1 Syngeneic Thy1.1 3` =  c('LLC1_StealTHY_Tum_3',	'LLC1_StealTHY_Met_3'),
   `LLC1 Syngeneic Thy1.1 4` =  c('LLC1_StealTHY_Tum_4',	'LLC1_StealTHY_Met_4'),
   `LLC1 Syngeneic Thy1.1 5` =  c('LLC1_StealTHY_Tum_5',	'LLC1_StealTHY_Met_5')
  ),
  `Comparison 3: CT26 Syngeneic Thy1.1` = list( 
    `CT26 Syngeneic Thy1.1 1` = c('ColoCT26_StealTHY_Tum_1',	'ColoCT26_StealTHY_Met_1'),
    `CT26 Syngeneic Thy1.1 2` = c('ColoCT26_StealTHY_Tum_2',	'ColoCT26_StealTHY_Met_2'),
    `CT26 Syngeneic Thy1.1 3` = c('ColoCT26_StealTHY_Tum_3',	'Replacement_1'),
    `CT26 Syngeneic Thy1.1 4` = c('Replacement_2',	'Replacement_3'),
    `CT26 Syngeneic Thy1.1 5` = c('Replacement_4',	'Replacement_4bis')
  ),
  `Comparison 4: CMT167 Syngeneic dCas9 + Puro` = list(
    `CMT167 Syngeneic dCas9 + Puro 3` = c('dC9BL6_LungIV_C167_3',	'dC9BL6_LiverIV_C167_3'),
    `CMT167 Syngeneic dCas9 + Puro 5` = c('dC9BL6_LungIV_C167_5',	'dC9BL6_LiverIV_C167_5')
  ),
  
  `Comparison 5: LLC1 Syngeneic dCas9 + Puro` = list(
    `LLC1 Syngeneic dCas9 + Puro 1` = c('dC9BL6_LungIV_LLC1_1',	'dC9BL6_LiverIV_LLC1_1'),
    `LLC1 Syngeneic dCas9 + Puro 3` = c('dC9BL6_LungIV_LLC1_3',	'dC9BL6_LiverIV_LLC1_3'),
    `LLC1 Syngeneic dCas9 + Puro 4` = c('dC9BL6_LungIV_LLC1_4',	'dC9BL6_LiverIV_LLC1_4'),
    `LLC1 Syngeneic dCas9 + Puro 5` = c('dC9BL6_LungIV_LLC1_5',	'dC9BL6_LiverIV_LLC1_5')
  ),
  
  `Comparison 6: CMT167 NSG dCas9 + Puro` = list(
    `CMT167 NSG dCas9 + Puro 1` = c('dC9NSG_LungIV_C167_1',	'dC9NSG_LiverIV_C167_1'),
    `CMT167 NSG dCas9 + Puro 2` = c('dC9NSG_LungIV_C167_2',	'dC9NSG_LiverIV_C167_2'),
    `CMT167 NSG dCas9 + Puro 3` = c('dC9NSG_LungIV_C167_3',	'dC9NSG_LiverIV_C167_3'),
    `CMT167 NSG dCas9 + Puro 4` = c('dC9NSG_LungIV_C167_4',	'dC9NSG_LiverIV_C167_4'),
    `CMT167 NSG dCas9 + Puro 5` = c('dC9NSG_LungIV_C167_5',	'dC9NSG_LiverIV_C167_5')
  ),
  
  `Comparison 7: LLC1 NSG dCas9 + Puro` = list(
    `LLC1 NSG dCas9 + Puro 1` = c('dC9NSG_LungIV_LLC1_1',	'dC9NSG_LiverIV_LLC1_1'),
    `LLC1 NSG dCas9 + Puro 2` = c('dC9NSG_LungIV_LLC1_2',	'dC9NSG_LiverIV_LLC1_2'),
    `LLC1 NSG dCas9 + Puro 3` = c('dC9NSG_LungIV_LLC1_3',	'dC9NSG_LiverIV_LLC1_3'),
    `LLC1 NSG dCas9 + Puro 4` = c('dC9NSG_LungIV_LLC1_4',	'dC9NSG_LiverIV_LLC1_4'),
    `LLC1 NSG dCas9 + Puro 5` = c('dC9NSG_LungIV_LLC1_5',	'dC9NSG_LiverIV_LLC1_5')
  ),
  `Comparison 10: Thy1.1 in BL6` = list(
    `CMT167 BL6 Thy1.1` = c('Replacement_9',	'Replacement_1bis'),
    `LLC1 BL6 Thy1.1` = c('Replacement_5',	'Replacement_6')
  ),
  `Comparison 11: Thy1.1 in NSG` = list(
    `CMT167 NSG Thy1.1` = c('Replacement_2bis',	'Replacement_3bis'),
    `LLC1 NSG Thy1.1` = c('Replacement_7',	'Replacement_8')
  )
)
```

Generate muller plots

``` r
cpm_threshold <- 10000
comp_list <- comp_list_muller

use_df <- assay(se, 'cpm') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('guide') %>% 
  mutate(
    gene = rowData(se)$Gene,
    gene = ifelse(grepl('Non_Target', gene), 'Non_Target', gene)
  ) %>%
  pivot_longer(-c(guide, gene), names_to = 'sample_alias', values_to = 'cpm')

j <-  names(comp_list)[1]
muller_plots <- foreach(j = names(comp_list)) %do% {
  i <- comp_list[[j]][[1]]
  res <- foreach(iname = names(comp_list[[j]])) %do% {
    i <- comp_list[[j]][[iname]]
    x <- use_df %>% 
      filter(sample_alias %in% i) %>% 
      mutate(
        sample_type = ifelse(sample_alias == i[1], 'Primary tumor', 'Metastasis'),
        sample_type = factor(sample_type, levels = c('Primary tumor', 'Metastasis'))
      )
    
    # Select guides with cpm > cpm_threshold
    x_guides_keep <- x %>% 
      group_by(sample_type, guide) %>% 
      summarise(cpm = max(cpm)) %>% 
      filter(cpm > cpm_threshold) %>% 
      pull(guide)
    
    # Group guides with cpm < cpm_threshold into the low abundant category, and sum counts for low abundant
    x %<>% 
      mutate(
        guide = ifelse(guide %in%  x_guides_keep, guide, 'low abundant')
      ) %>% 
      group_by(sample_type, guide) %>%
      summarise(cpm = sum(cpm)) %>% 
      ungroup() %>%
      mutate(
        guide = fct_reorder(guide, cpm),
        guide = relevel(guide, ref = 'low abundant')
      )
    
    use_cols <- c(
      'low abundant' = 'grey80',
      colorRampPalette(rev(palette_OkabeIto))(nlevels(x$guide) - 1) %>% set_names(levels(x$guide)[-1])
    )
    
    x %>% 
      ggplot( aes(x = sample_type, y = cpm/10000, group = guide, fill = guide)) + 
      geom_area(colour = alpha("white", 0.1), linewidth = 0.08, alpha = 0.8) +
      scale_fill_manual(values = use_cols, guide = guide_legend(ncol = 3)) +
      labs(
        x = '',
        y = expression("CPM x 10"^-4),
        title = iname,
        fill = ''
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1000000/10000)) +
      guides(fill = guide_legend(nrow=8)) +
      theme(
        panel.background = element_rect(fill = "grey90"),
        # plot.margin = margin(0.5, 1.5, 0.5, 0, "cm"),
        plot.margin = margin(0, 0.5, 0, 0, "cm"),
        axis.text =  element_text(size=3),
        axis.title =  element_text(size=3),
        plot.title = element_text(size=3, hjust = 0.5),
        legend.text = element_text(size=3),
        legend.title = element_text(size=3),
        legend.position = 'bottom',
        legend.key.size = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0.1, 'mm'),
        legend.key.spacing.x = unit(0.1, 'mm')
      )
  }
  names(res) <- names(comp_list[[j]])
  return(res)

}
names(muller_plots) <- names(comp_list)
# muller_syngeneic_representative <- muller_plots
```

### Supplementary Figure 2D (left upper)

Muller plots depicting the relative frequency of unique sgRNAs in matched primary tumors (mammary or lung orthotopic) versus synchronous metastasis samples (lung metastasis for mammary tumors, liver metastasis for lung tumors) from syngeneic and NSG recipients transplanted with the indicated models previously transduced as indicated; n=5 per group (representative samples are shown).


``` r
plot_list <- list(
  muller_plots[['Comparison 7: LLC1 NSG dCas9 + Puro']][['LLC1 NSG dCas9 + Puro 2']] + guides(fill = 'none'),
  muller_plots[['Comparison 5: LLC1 Syngeneic dCas9 + Puro']][['LLC1 Syngeneic dCas9 + Puro 1']]+  guides(fill = 'none'),
  muller_plots[['Comparison 6: CMT167 NSG dCas9 + Puro']][['CMT167 NSG dCas9 + Puro 4']] + guides(fill = 'none') ,
  muller_plots[['Comparison 4: CMT167 Syngeneic dCas9 + Puro']][['CMT167 Syngeneic dCas9 + Puro 3']] + guides(fill = 'none')
)

wrap_plots(plot_list, ncol = 2, byrow = TRUE) %>% print
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig2-d-left-upper-1.png" style="display: block; margin: auto;" />

### Supplementary Figure 2D (right upper)

Muller plots depicting the relative frequency of unique sgRNAs in matched primary tumors (mammary or lung orthotopic) versus synchronous metastasis samples (lung metastasis for mammary tumors, liver metastasis for lung tumors) from syngeneic and NSG recipients transplanted with the indicated models previously transduced as indicated; n=5 per group (representative samples are shown).


``` r
plot_list <- list(
  muller_plots[['Comparison 11: Thy1.1 in NSG']][['LLC1 NSG Thy1.1']]+ guides(fill = 'none'),
  muller_plots[['Comparison 10: Thy1.1 in BL6']][['LLC1 BL6 Thy1.1']]+ guides(fill = 'none'),
  muller_plots[['Comparison 11: Thy1.1 in NSG']][['CMT167 NSG Thy1.1']]+ guides(fill = 'none'),
  muller_plots[['Comparison 10: Thy1.1 in BL6']][['CMT167 BL6 Thy1.1']]+ guides(fill = 'none')
)

wrap_plots(plot_list, ncol = 2, byrow = TRUE) %>% print
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig2-d-right-upper-1.png" style="display: block; margin: auto;" />

### Supplementary Figure 2F

Muller plot depicting the relative frequency of the predominant sgRNAs in matched primary tumors versus synchronous metastasis samples (lungs) from syngeneic recipients transplanted with CT26 cells, previously transduced with Thy1 as indicated; n=5 each (representative samples are shown).


``` r
print(muller_plots[['Comparison 3: CT26 Syngeneic Thy1.1']][['CT26 Syngeneic Thy1.1 3']] + guides(fill = 'none'))
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig2-f-1.png" style="display: block; margin: auto;" />

### Supplementary Figure 4L

Application of StealTHY to syngeneic models of lung and colon carcinoma. Left: schematic representation of the in vivo CRISPR screen strategy to identify metastatic genes in immunocompetent mice, showing the experimental pipeline followed for each cancer type. Right: Muller plots depicting the relative frequency of the predominant sgRNAs in matched primary tumors (lung or colon orthotopic) versus synchronous metastasis samples (liver metastasis for lung tumors, lung metastasis for colon tumors) from syngeneic recipients transplanted with the indicated cell models, which received StealTHY KO of the mouse interactome library; n=5 each (representative samples are shown).


``` r
plot_list <- list(
  muller_plots[['Comparison 2: LLC1 Syngeneic Thy1.1']][['LLC1 Syngeneic Thy1.1 4']],
  muller_plots[['Comparison 3: CT26 Syngeneic Thy1.1']][['CT26 Syngeneic Thy1.1 1']]
)

wrap_plots(plot_list, ncol = 1) %>% print
```

<img src="figure/crispr-mm_2215_sgRNA-clonality.Rmd/muller-sfig4-l-1.png" style="display: block; margin: auto;" />





<br>
<p>
<button type="button" class="btn btn-default btn-workflowr btn-workflowr-sessioninfo"
  data-toggle="collapse" data-target="#workflowr-sessioninfo"
  style = "display: block;">
  <span class="glyphicon glyphicon-wrench" aria-hidden="true"></span>
  Session information
</button>
</p>

<div id="workflowr-sessioninfo" class="collapse">

``` r
sessionInfo()
```

R version 4.4.3 (2025-02-28)
Platform: aarch64-apple-darwin20
Running under: macOS 26.0.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Zurich
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] cowplot_1.2.0               ggpubr_0.6.1               
 [3] ggbeeswarm_0.7.2            colorblindr_0.1.0          
 [5] colorspace_2.1-1            patchwork_1.3.1            
 [7] ggh4x_0.3.1                 ggridges_0.5.6             
 [9] SummarizedExperiment_1.36.0 Biobase_2.66.0             
[11] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[13] IRanges_2.40.1              S4Vectors_0.44.0           
[15] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
[17] matrixStats_1.5.0           diptest_0.77-1             
[19] kableExtra_1.4.0            DT_0.34.0                  
[21] magrittr_2.0.3              foreach_1.5.2              
[23] knitr_1.50                  lubridate_1.9.4            
[25] forcats_1.0.0               stringr_1.5.1              
[27] dplyr_1.1.4                 purrr_1.1.0                
[29] readr_2.1.5                 tidyr_1.3.1                
[31] tibble_3.3.0                ggplot2_3.5.2              
[33] tidyverse_2.0.0             workflowr_1.7.1            

loaded via a namespace (and not attached):
 [1] rlang_1.1.6             git2r_0.36.2            compiler_4.4.3         
 [4] getPass_0.2-4           systemfonts_1.2.3       callr_3.7.6            
 [7] vctrs_0.6.5             pkgconfig_2.0.3         crayon_1.5.3           
[10] fastmap_1.2.0           backports_1.5.0         XVector_0.46.0         
[13] labeling_0.4.3          promises_1.3.3          rmarkdown_2.29         
[16] tzdb_0.5.0              UCSC.utils_1.2.0        ps_1.9.1               
[19] xfun_0.53               zlibbioc_1.52.0         cachem_1.1.0           
[22] jsonlite_2.0.0          later_1.4.4             DelayedArray_0.32.0    
[25] broom_1.0.9             R6_2.6.1                bslib_0.9.0            
[28] stringi_1.8.7           RColorBrewer_1.1-3      car_3.1-3              
[31] jquerylib_0.1.4         Rcpp_1.1.0              iterators_1.0.14       
[34] httpuv_1.6.16           Matrix_1.7-3            timechange_0.3.0       
[37] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1      
[40] abind_1.4-8             yaml_2.3.10             codetools_0.2-20       
[43] processx_3.8.6          lattice_0.22-7          withr_3.0.2            
[46] evaluate_1.0.5          xml2_1.3.8              pillar_1.11.0          
[49] carData_3.0-5           whisker_0.4.1           generics_0.1.4         
[52] rprojroot_2.1.0         hms_1.1.3               scales_1.4.0           
[55] glue_1.8.0              tools_4.4.3             ggsignif_0.6.4         
[58] fs_1.6.6                grid_4.4.3              crosstalk_1.2.2        
[61] GenomeInfoDbData_1.2.13 beeswarm_0.4.0          vipor_0.4.7            
[64] Formula_1.2-5           cli_3.6.5               textshaping_1.0.1      
[67] S4Arrays_1.6.0          viridisLite_0.4.2       svglite_2.2.1          
[70] gtable_0.3.6            rstatix_0.7.2           sass_0.4.10            
[73] digest_0.6.37           SparseArray_1.6.2       htmlwidgets_1.6.4      
[76] farver_2.1.2            htmltools_0.5.8.1       lifecycle_1.0.4        
[79] httr_1.4.7             
</div>
