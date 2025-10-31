---
title: "Analysis of AMHR2 KO RNA-seq data"
author: "Francesc Castro-Giner"
date: 'October 31, 2025'
output: workflowr::wflow_html
editor_options:
  chunk_output_type: console
params:
  date: 'October 31, 2025'
  se_path: ./output/rnaseq/amhr2_ko/se.rds
  output_dir: ./output/rnaseq/amhr2_ko
  min_counts: 5
  min_present_prop: 0.40
  msigdb_dir: ./data/resources/MSigDB
  msigdb_vs: v2023.2.Mm
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
	Untracked:  analysis/crispr-hsapiens_2180_sgRNA.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-StealTHY.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-clonality.md
	Untracked:  analysis/index.md
	Untracked:  analysis/rnaseq-tumor-bulk.md
	Untracked:  analysis/style.css
	Untracked:  analysis/tcga_survival_crispr_hits_signatures.md
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
to the R Markdown (<code>analysis/rnaseq-amhr2_ko-deg.Rmd</code>) and HTML (<code>docs/rnaseq-amhr2_ko-deg.html</code>)
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
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/455fa9121fcca51363cf639f7ee170fb108e34c1/analysis/rnaseq-amhr2_ko-deg.Rmd" target="_blank">455fa91</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Prepare files for publication</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/8632845b1c6b024330836a1a798537363e31980e/docs/rnaseq-amhr2_ko-deg.html" target="_blank">8632845</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/bc402e6695932824d0587de6c03eb192be28dfca/analysis/rnaseq-amhr2_ko-deg.Rmd" target="_blank">bc402e6</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>uodate rnaseq analysis</td>
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
knitr::opts_chunk$set(results='asis', echo=TRUE, message=FALSE, warning=FALSE, error=FALSE, fig.align = 'center', fig.width = 3.5, fig.asp = 0.618, dpi = 600, dev = c("png", "pdf"), engine.opts = list(bash = "-l"))

options(stringsAsFactors = FALSE)

use_seed <- 1100101
set.seed(use_seed)

dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
```

Load packages

``` r
library(tidyverse)
library(showtext)
library(foreach)
library(DT)
library(knitr)
library(kableExtra)
library(cowplot)
library(colorblindr)
library(RColorBrewer)
library(ggbeeswarm)
library(ggrepel)
library(magrittr)
library(openxlsx)

library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(scuttle)
library(scran)
library(pheatmap)
library(PCAtools)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(circlize)
```


Load ggplot theme

``` r
source("./configuration/rmarkdown/ggplot_theme.R")
```

Set color palettes


Load custom functions

``` r
source('./code/R-functions/dge_wrappers.r')
source('./code/R-functions/dge_report.r')
source('./code/R-functions/gse_omnibus.r')
source('./code/R-functions/gse_report.r')
source('./code/R-functions/subchunkify.R')
```

Load MSigDB gene sets

``` r
gmt_files_symbols <- list(
  msigdb.h = paste0(params$msigdb_dir, '/', params$msigdb_vs, '/mh.all.', params$msigdb_vs, '.symbols.gmt')
)
```

Clean files generated in previous runs

``` r
rmd_file <- current_input()
if(!is.null(rmd_file)) {
  figures_dir <- file.path('./docs/figure',rmd_file)
  if(dir.exists(figures_dir)) {
    unlink(file.path(figures_dir, "*"))
  }
}
```


## Configure analyses

### Data wrangling
Load SE object

``` r
se <- readRDS(params$se_path)
```

Modify sample annotation

``` r
colData(se) %<>% 
  data.frame %>% 
  mutate(
    condition = sample_alias %>% gsub('_rep.', '', .),
    condition_short = case_match(
      condition,
      "AMHR2_Control_culture" ~ "c_AMHR2_ctrl",
      "AMHR2_KO_culture" ~ "c_AMHR2_KO",
      "AMHR2_KO_noTGFbeta" ~ "AMHR2_KO_noTGFb",
      "AMHR2_Overexpression" ~ "AMHR2_over",
      "InVivo_AMHR2_KO" ~ "i_AMHR2_KO",
      "InVivo_Control_sgRNA" ~ "i_ctrl_sgRNA"
    ),
    sample_type = condition_short,
    condition = factor(condition),
    condition_short = factor(condition_short)

  ) %>% DataFrame
```


Filtering out low-abundance features. Keeping those features with at least 10 counts in at least 50% of samples of the smallest group.

``` r
use_rows <- filterByExpr(se, group=se$condition, min.count = 10, min.total.count = 15, min.prop = 0.5)
se <- se[use_rows,]
```

TMM normalization

``` r
norm_method <- "TMM"
exprs_mat_for_norm <- assay(se, 'counts')
norm_factors <- calcNormFactors(exprs_mat_for_norm, method = norm_method)
lib_size <- colSums(exprs_mat_for_norm)/10^6
se$size_factor <- (norm_factors * lib_size)
se <- logNormCounts(se, size.factors =se$size_factor)
se <- logNormCounts(se, size.factors =se$size_factor, log = FALSE)
```

CPM calculation

``` r
assay(se, 'cpm') <- edgeR::cpm(assay(se, 'normcounts'))
assay(se, 'logcpm') <- edgeR::cpm(assay(se, 'normcounts'), log = TRUE)
```


### Define of comparisons

#### Process data

``` r
x <- colData(se) %>% data.frame

# Define groups
use_levels <- levels(se$condition_short)
i <- use_levels[1]
group_samples <- foreach(i = use_levels) %do% {
  use_cols <- se$condition_short == i
  colnames(se)[use_cols]
}
names(group_samples) <- use_levels


# Define all potential pairwise comparisons
group_combn <- combn(names(group_samples) , 2, simplify = FALSE)

i <- group_combn[[1]]
comp_list <- foreach(i = group_combn) %do%{
  res <- list(
    group_samples[[i[1]]],
    group_samples[[i[2]]]
  )
  names(res) <- i
  res
}
names(comp_list) <- lapply(group_combn, function(x) paste0(x[1], '-o-', x[2]))


# Add additional comparison
comp_list[['AMHR2_KO-o-ctrl']] <- list(
  `AMHR2_KO` = x %>% 
    filter(condition %in% c('InVivo_AMHR2_KO', 'AMHR2_KO_culture')) %>% 
     pull(sample_alias),
  `control` = x %>% 
    filter(condition %in% c('InVivo_Control_sgRNA', 'AMHR2_Control_culture')) %>% 
     pull(sample_alias)
)

keep_comp <- c(
  `InVivo_AMHR2_KO + AMHR2_KO_culture--over--InVivo_Control_sgRNA + AMHR2_Control_culture` = 'AMHR2_KO-o-ctrl',
  `AMHR2_KO_noTGFbeta--over--AMHR2_Control_culture` = 'AMHR2_KO_noTGFb-o-c_AMHR2_ctrl',
  `InVivo_AMHR2_KO–over–InVivo_Control_sgRNA` = 'i_AMHR2_KO-o-i_ctrl_sgRNA'
)
comp_list <- comp_list[keep_comp]


# Generate DF with name description
comp_list_description <- data.frame(
  comparison = names(comp_list),
  comparison_long = names(keep_comp)
)
```

#### Table of comparisons

We filter comparisons with n < 3 samples in each group

``` r
i = names(comp_list)[1]
comp_df <- foreach(i = names(comp_list), .combine = rbind) %do% {
  icomp <- comp_list[[i]]
  data.frame(
    comparison = i,
    `Sample Group (n)` =  length(icomp[[1]]),
    `Reference Group (n)` = length(icomp[[2]]),
    check.names = FALSE
  )
}

# filter comparisons were we have at least 3 samples in each group
comp_df %<>% 
  mutate(Evaluated = `Sample Group (n)` >= 3 & `Reference Group (n)` >= 3)
kepp_comp <- comp_df %>% filter(Evaluated) %>% pull(comparison)
comp_list <- comp_list[kepp_comp]

# Add long names
comp_df <- comp_df %>% left_join(comp_list_description) %>% 
  mutate(comparison = comparison_long) %>%
  dplyr::select(-comparison_long)

comp_df %>%
  datatable(., 
            rownames = FALSE, 
            filter = 'top', 
            caption = 'List of comparisons for differential expression. Comparisons with n < 3 samples in each group won\'t be evaluated',
            extensions = 'Buttons', 
            options = list(
              dom = 'Blfrtip',
              buttons = c('csv', 'excel')
            ))
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-0e41329c1963262c826a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0e41329c1963262c826a">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"3\" data-max=\"6\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"3\" data-max=\"6\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"logical\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"width: 100%; display: none;\">\n      <select multiple=\"multiple\" style=\"width: 100%;\" data-options=\"[&quot;true&quot;,&quot;false&quot;]\"><\/select>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>List of comparisons for differential expression. Comparisons with n &lt; 3 samples in each group won't be evaluated<\/caption>","data":[["InVivo_AMHR2_KO + AMHR2_KO_culture--over--InVivo_Control_sgRNA + AMHR2_Control_culture","AMHR2_KO_noTGFbeta--over--AMHR2_Control_culture","InVivo_AMHR2_KO–over–InVivo_Control_sgRNA"],[6,3,3],[6,3,3],[true,true,true]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>comparison<\/th>\n      <th>Sample Group (n)<\/th>\n      <th>Reference Group (n)<\/th>\n      <th>Evaluated<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"className":"dt-right","targets":[1,2]},{"name":"comparison","targets":0},{"name":"Sample Group (n)","targets":1},{"name":"Reference Group (n)","targets":2},{"name":"Evaluated","targets":3}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>
```


## Data generation

### Run differential expression analysis
Get Ensembl gene info

``` r
ensembl <-  biomaRt::useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl",
                                mirror = 'www')

gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'),
                              filters = 'external_gene_name',
                              values = rowData(se)$gene_name, mart =ensembl) %>%
    dplyr::rename('gene_name' = 'external_gene_name') %>%
    unique


gene_desc <- gene_desc[!duplicated(gene_desc$gene_name),]
saveRDS(gene_desc, file = file.path(dirname(params$se_path), 'ensembl_gene_desc.rds'))
```


``` r
gene_desc <- readRDS(file = file.path(dirname(params$se_path), 'ensembl_gene_desc.rds'))
```

Run differential expression analysis using EdgeR QLF

``` r
use_comp <-  comp_list[[1]]
dge_list <- foreach(use_comp = comp_list) %do%{
  se_cols <- se$sample_alias %in% unlist(use_comp)
  use_se <- se[,se_cols]
  use_se$group <- ifelse(
    use_se$sample_alias %in% use_comp[[1]],
    names(use_comp)[1],
    names(use_comp)[2]
  )
  
  if(is.null(use_comp$batch)) {
    dge <- edgeR_dge(
      use_se,
      # Design configuration for differential expression
      group_var =  'group',
      group_sample = names(use_comp)[1],
      group_ref = names(use_comp)[2],
      design_formula = "~ group",
      coef = 'last',
      # Conversion from SingleCellExperiment/Summarized Experiment to DGEList
      assay_to_DGEList = 'counts',
      assay_to_row_filter = "counts",
      # Feature filtering parameters
      use_filterByExpr = TRUE,
      min_counts = params$min_counts,
      min_present_prop = params$min_present_prop,
      # EdgeR workflow configuration
      run_calcNormFactors = 'TMM',
      estimateDisp_robust = TRUE,
      estimateDisp_trend.method = "locfit",
      glm_approach = "QLF",
      # Output configuration
      adjust_method = 'BH'
      )
    } else {
      design_formula = paste0("~ ", use_comp$batch," + group")
      dge <- edgeR_dge(
        use_se,
        # Design configuration for differential expression
        group_var =  'group',
        group_sample = names(use_comp)[1],
        group_ref = names(use_comp)[2],
        batch_vars = use_comp$batch,
        design_formula = design_formula,
        coef = 'last',
        # Conversion from SingleCellExperiment/Summarized Experiment to DGEList
        assay_to_DGEList = 'counts',
        assay_to_row_filter = "counts",
        # Feature filtering parameters
        use_filterByExpr = TRUE,
        min_counts = params$min_counts,
        min_present_prop = params$min_present_prop,
        # EdgeR workflow configuration
        run_calcNormFactors = 'TMM',
        estimateDisp_robust = TRUE,
        estimateDisp_trend.method = "locfit",
        glm_approach = "QLF",
        # Output configuration
        adjust_method = 'BH'
        )
    }
  
  # Add gene description
  use_res <- dge$results %>%  left_join(., gene_desc, relationship = "many-to-many")

  dge$results <- use_res %>% 
    filter(!duplicated(feature)) %>% 
    mutate(rownames = feature) %>% 
    column_to_rownames('rownames')
  
   return(dge)
}

names(dge_list) <- names(comp_list)

# Save objects
saveRDS(dge_list, file = file.path(params$output_dir, 'dge-edgeR_QLF.rds'))
```

### Run functional enrichment

``` r
dge_list <- readRDS(file.path(params$output_dir, 'dge-edgeR_QLF.rds'))
dge <-  dge_list[[2]]

gse_list <- foreach(dge = dge_list) %do%{
  gse_res <- gse_omnibus(
    feature_names = dge$results$gene_name,
    p = dge$results$FDR,
    fc = dge$results$logFC,
    fc_thrs = 0.5,
    gmt_files = gmt_files_symbols,
    
    annot_db = 'org.Mm.eg.db',
    organism = 'mmusculus',
    kegg_organism = 'mmu',
    reactome_organism = 'mouse',
    
    run_enricher = FALSE,
    run_enrichGO = FALSE,
    run_enrichKEGG = FALSE,
    run_enrichReactome = FALSE,
    
    run_GSEA = TRUE,
    run_gseGO = FALSE,
    run_gseKEGG = FALSE,
    run_gseReactome = FALSE,
    
    go_ontologies = c('BP', 'MF'),
    go_simplify_p = 0.05,
    args_gse = list(minGSSize = 10, maxGSSize = 250, 
                    pvalueCutoff = 1, 
                    BPPARAM = BiocParallel::MulticoreParam(1))
    )
  
   return(gse_res)

}
names(gse_list) <- names(dge_list)

# Save objects
saveRDS(gse_list, file = file.path(params$output_dir,'gse_omnibus-edgeR_QLF.rds'))
```


### Load data

``` r
dge_list <- readRDS(file = file.path(params$output_dir, 'dge-edgeR_QLF.rds'))
gse_list <- readRDS(file = file.path(params$output_dir, 'gse_omnibus-edgeR_QLF.rds'))
```


### Knitr:template

``` r
selected_labels_volcano <- c('Amhr2')
# knit template: https://workflowr.github.io/workflowr/articles/wflow-07-common-code.html
src <- knit_expand(
  file = here::here("analysis/templates/rnaseq-deg-edger.Rmd")
  )
```
<!-- To be compatible with workflowr, we need to specify the absolute directory under which the plots are generated, otherwise they will be saved into figure/<basename>.Rmd/<chunk-name>.<extension. and not visible in the html document: https://github.com/workflowr/workflowr/issues/111 -->



## Download results

<!-- Write DGE results to xlsx file -->

``` r
# File name summary
rmd_file <- current_input()
if(is.null(rmd_file))
  rmd_file <- 'tmp'
file_xlsx <- file.path('./docs/file',rmd_file, 'dge.xlsx')

dir.create(dirname(file_xlsx), recursive = TRUE, showWarnings = FALSE)

# Generate workbook
wb <- createWorkbook()
for(i in names(dge_list)) {
  addWorksheet(wb, i)
  res <- dge_list[[i]]$results %>% 
    dplyr::select(feature, gene_name, mean, detected, logFC, PValue, FDR)
  writeData(wb, i, res)
}
saveWorkbook(wb, file_xlsx, TRUE)

dge_file_xlsx <- file_xlsx
```

<!-- Write GSE results to xlsx file -->

``` r
gse_collections <- list(
  GSEA = names(gmt_files_symbols),
  gseGO = c("BP", "MF")
)

# File name summary
rmd_file <- current_input()
if(is.null(rmd_file))
  rmd_file <- 'tmp'

file_xlsx_list <- list()
for(sr in names(gse_collections)) {
  file_xlsx_list[[sr]] <- list()
  for(gc in gse_collections[[sr]]) {
    # Generate file
    file_xlsx <- file.path('./docs/file', 
                           rmd_file, 
                           paste0('gse_', tolower(sr), '_', tolower(gc), '.xlsx')
                           )
    file_xlsx_list[[sr]][[gc]] <- file_xlsx
    dir.create(dirname(file_xlsx), recursive = TRUE, showWarnings = FALSE)
    
    # Generate workbook
    wb <- createWorkbook()
    for(i in names(gse_list)) {
      addWorksheet(wb, i)
      if(!is.null(gse_list[[i]][[sr]][[gc]])) {
        res <- gse_list[[i]][[sr]][[gc]]@result
        writeData(wb, i, res)
      }
    }
    saveWorkbook(wb, file_xlsx, TRUE)
  }
}
```

<!-- Change dge_list names after writing to xlsx -->

``` r
comp_names <- comp_list_description$comparison_long %>% 
  set_names(comp_list_description$comparison)
names(comp_list) <- comp_names[names(comp_list)]
names(dge_list) <- comp_names[names(dge_list)]
names(gse_list) <- comp_names[names(gse_list)]
```


The tables of results can be downloaded using the following links:

- [**Differential gene expression**](./file/rnaseq-amhr2_ko-deg.Rmd/dge.xlsx)
<!-- Write links to download GSE results -->
- [**Gene set enrichment analysis: msigdb.h**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsea_msigdb.h.xlsx) 

- [**Gene set enrichment analysis: GO:BP**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsego_bp.xlsx) 

- [**Gene set enrichment analysis: GO:MF**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsego_mf.xlsx) 

## Functional enrichment analysis

``` r
gse_collections <- list(
  GSEA = names(gmt_files_symbols)
)
```

### Download results
<!-- Write links to download GSE results -->
- [**Gene set enrichment analysis: msigdb.h**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsea_msigdb.h.xlsx) 

- [**Gene set enrichment analysis: GO:BP**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsego_bp.xlsx) 

- [**Gene set enrichment analysis: GO:MF**](./file/rnaseq-amhr2_ko-deg.Rmd/gse_gsego_mf.xlsx) 

### Summary of results

``` r
res <- foreach(i = names(gse_list), .combine = rbind) %do% {
  foreach(sr = names(gse_collections), .combine = rbind) %do% {
    foreach(gc = gse_collections[[sr]], .combine = rbind) %do% {
      if(sr == 'GSEA') {
        x <- gse_list[[i]][[sr]][[gc]]@result %>%
          filter(p.adjust <= 0.05)
        if(gc == 'msigdb.c2.cp') {
          x %<>%
            filter(grepl('REACTOME|KEGG|WP', ID))
        }
      }
      if(sr == 'gseGO') {
        x <- gse_list[[i]][[sr]][[gc]]@result %>%
          filter(p.adjust <= 0.05)
      }
      c(
        Comparison = i,
        `Collection` = gc,
        `Total` = nrow(x),
        `Supressed` = x %>% filter(NES < 0) %>% nrow,
        `Activated` = x %>% filter(NES > 0) %>% nrow
      )
    }
  }
}
res <- res %>% data.frame(check.names = F)
```

### Dotplots
<!-- # Subchunkify dotplots -->
<!-- ## chunk options need to be set to echo=FALSE, results='asis' -->
<!-- ## Cause results='asis',we must use raw HTML to print text and headings -->
<h4> <span class="header-section-number"> 5.3.1 </span>Comparison InVivo_AMHR2_KO + AMHR2_KO_culture--over--InVivo_Control_sgRNA + AMHR2_Control_culture </h4>

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1">
  Past versions of gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1.png
  </button>
  </p>

  <div id="fig-gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/8632845b1c6b024330836a1a798537363e31980e/docs/figure/rnaseq-amhr2_ko-deg.Rmd/gse-dotplot-invivo_amhr2_ko-+-amhr2_ko_culture--over--invivo_control_sgrna-+-amhr2_control_culture_GSEA_msigdb.h-1.png" target="_blank">8632845</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  
                     

<h4> <span class="header-section-number"> 5.3.2 </span>Comparison AMHR2_KO_noTGFbeta--over--AMHR2_Control_culture </h4>

<h4> <span class="header-section-number"> 5.3.3 </span>Comparison InVivo_AMHR2_KO–over–InVivo_Control_sgRNA </h4>

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gse-dotplot-invivo_amhr2_ko–over–invivo_control_sgrna_GSEA_msigdb.h-1.png" style="display: block; margin: auto;" />
                     









## GSEA plots for selected signatures

Define signatures

``` r
manual_gene_sets <- list(
   # Signature 4
  `STAT1 target genes and interactors` = c(
    'Ciita',
    'Socs1',
    'Cxcl10',
    'Irf1',
    'Irf7',
    'Irf9',
    'Stat1',
    'Stat2',
    'Stat3',
    'Stat5a',
    'Stat5b',
    'Sting1',
    'Jak1',
    'Jak2',
    'Tap1',
    'Tap2',
    'Usp18',
    'Cyp1b1',
    'Mx1',
    'Ido1',
    'Bcl6',
    'Gbp5',
    'Egr1',
    'Fos',
    'Lyn'
  ),
  # Signature 2
  `EMT regulators and targets` = c(
    'Epcam',
    'Cdh1',
    'Cdh2',
    'Cldn1',
    'Cldn3',
    'Cldn4',
    'Cldn7',
    'Vim',
    'Acta2',
    'Snai2',
    'Twist1',
    'Zeb1',
    'Zeb2',
    'Grhl2',
    'Esrp1',
    'Esrp2',
    'Marveld2',
    'Rab25',
    'Myo5C',
    'Pkp1',
    'Pkp3',
    'Jup',
    'Dsp',
    'Ppl',
    'Dsg3'
  ),
  # Signature 1
  `Mammary epithelium proliferation` = c(
    'Aurka',
    'Aurkb',
    'Brca1',
    'Brca2',
    'Ccdn1',
    'Cdkn2a',
    'Esr1',
    'Gata3',
    'Hoxa5',
    'Mki67',
    'Plk1',
    'Wnt5A',
    'Agap2',
    'Areg',
    'Bax',
    'Deaf1',
    'Epha2',
    'Etv4',
    'Id2',
    'Kdm5b',
    'Robo1',
    'Mapk1',
    'Med1',
    'Pygo2',
    'Tfap2c'
  ),
  # Signature 3
  `SMAD target genes and interactors` = c(
    'Wnt4',
    'Cited2',
    'Tob1',
    'Acvr1',
    'Hipk2',
    'Ppm1a',
    'Crebbp',
    'Zfyve9',
    'Col5a2',
    'Tgif1',
    'Ankrd1',
    'Smad4',
    'Smurf1',
    'Magi2',
    'Zmiz1',
    'Trim33',
    'Skil',
    'Ipo7',
    'Mef2a',
    'Smad9',
    'Usp9x',
    'Rnf111',
    'Tcf12',
    'Hmga2',
    'Foxa2'
  )
)


i <- names(manual_gene_sets)[1]
manual_term2gene <- foreach(i = names(manual_gene_sets), .combine = rbind) %do% {
  cbind(
    term = rep(i, length(manual_gene_sets[[i]])), 
    gene = manual_gene_sets[[i]]
    )
}
```

Run functional enrichment

``` r
library(clusterProfiler)
dge <-  dge_list[[2]]

gse_manual_list <- foreach(dge = dge_list) %do%{
  feature_fc_list <- dge$results$logFC %>% 
    set_names(make.names(dge$results$gene_name, unique = TRUE)) %>% 
    sort(decreasing = TRUE)
  res <- GSEA(feature_fc_list, 
              TERM2GENE=manual_term2gene,
              minGSSize = 0,
              maxGSSize = Inf,
              pvalueCutoff = 1, verbose = FALSE)
  res
}

names(gse_manual_list) <- names(dge_list)

# Save objects
saveRDS(gse_manual_list, file = file.path(params$output_dir,'gsea-edgeR_QLF-signatures.rds'))
```


``` r
gse_manual_list <- readRDS(file = file.path(params$output_dir,'gsea-edgeR_QLF-signatures.rds'))
comp_names <- comp_list_description$comparison_long %>% 
  set_names(comp_list_description$comparison)
names(gse_manual_list) <- comp_names[names(gse_manual_list)]
```


``` r
gseaScores <- getFromNamespace("gseaScores", "DOSE")
gsInfo <- getFromNamespace("gsInfo", "enrichplot")
gseaplot2_custom <- function (
    x, 
    geneSetID, 
    title = "", 
    color = "green", 
    base_size = 11, 
    rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
    ES_geom = "line",
    axis.linewidth = (1/2.141959)/4,
    pres.linewidth = 1,
    p2.fill = 'transparent',
    p2.linerange.linewidth = 0.5,
    legend.position = c(0.8, 0.8),
    p.pos.ylab = 'Ranked List Metric',
    p2.labels = FALSE, # Only works when there are two geneSetID
    reverse.rank = FALSE,
    p2.add.intervals.colors = FALSE,
    p2.segment.color.var = 'Description', # Either Description or geneList (log2FC)
    geneList.limit = NULL,
    fontface = "italic" 
    )
{
    ##
    ## Enrichmment score plot
    ##
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    } else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    
    # Reverse rank
    if(reverse.rank) {
      gsdata$x <- rev(gsdata$x)
    }
    
    # Cap the geneList values
    if(!is.null(geneList.limit)) {
      gsdata$geneListCapped <- gsdata$geneList
      gsdata$geneListCapped[gsdata$geneListCapped > geneList.limit] <- geneList.limit
      gsdata$geneListCapped[gsdata$geneListCapped < -geneList.limit] <- -geneList.limit
      p2.segment.color.var <- 'geneListCapped'
    }
    
    p <- ggplot(gsdata, aes_(x = ~x)) + 
      xlab(NULL) + 
      theme_classic(base_size) + 
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()
      ) + 
      scale_x_continuous(expand = c(0, 0))
    
    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
            linewidth = 1)
    } else {
        es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
            linewidth = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p + 
      es_layer + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      theme(
        legend.position = legend.position,
        legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent")
        ) + 
      ylab("Running Enrichment Score") + 
      theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.border = element_rect(fill ="transparent",color = "black", linewidth = axis.linewidth)
        )
    ##
    ## Plot gene position in rank
    ##
    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == 
            term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    
    p2 <- ggplot(gsdata, aes_(x = ~x)) + 
      annotate("rect", xmin=0, xmax=100+max(gsdata$x), ymin=0 , ymax=2, alpha=0.8, color=p2.fill, fill=p2.fill) +
      geom_linerange(aes_string(ymin = 'ymin', ymax = 'ymax', color = p2.segment.color.var), linewidth = p2.linerange.linewidth) + 
      xlab(NULL) + 
      ylab(NULL) + 
      theme_classic(base_size) + 
      theme(
        legend.position = "none", 
        plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(fill ="transparent",color = "black", linewidth = axis.linewidth)
        ) + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    
    if(p2.labels == TRUE & length(geneSetID) == 2) {
      gsdata_labels <- gsdata %>% filter(ymax > 0)
      max_y <- max(gsdata_labels$ymax)
      
      description_labels <- unique(gsdata_labels$Description) 
      description_cut_point <- ceiling(length(description_labels) / 2)
      description_bottom <- description_labels[1:description_cut_point]
      description_top <- description_labels[(description_cut_point+1):length(description_labels)]
      
      gsdata_labels_bottom <- gsdata_labels %>% filter(Description %in% description_bottom)
      gsdata_labels_top <- gsdata_labels %>% filter(Description %in% description_top)
      max_y <- max(gsdata_labels$ymax)
      p2 <- p2 + 
        geom_text_repel(data = gsdata_labels_top,
                        aes(y = max_y, label = gene),
                        force_pull   = 0, # do not pull toward data points
                        nudge_y      = max_y * 0.15,
                        direction    = "x",
                        ylim = c(2.5, 4),
                        angle        = 90,
                        hjust        = 0,
                        max.time = 3,
                        max.overlaps = Inf,
                        min.segment.length = 0,
                        segment.size = 0.2,
                        segment.square = TRUE,
                        segment.inflect = FALSE,
                        segment.curvature = 1,
                        segment.angle = -1,
                        box.padding = 0.075,
                        fontface = "italic"
        ) +
        geom_text_repel(data = gsdata_labels_bottom,
                        aes(y = 0, label = gene),
                        force_pull   = 0, 
                        nudge_y      = -(max_y * 0.15),
                        direction    = "x",
                        ylim = c(-2, -0.5),
                        angle        = 90,
                        hjust        = 1,
                        max.time = 3,
                        max.overlaps = Inf,
                        min.segment.length = 0,
                        segment.size = 0.2,
                        segment.square = TRUE,
                        segment.inflect = FALSE,
                        segment.curvature = 1,
                        segment.angle = -1,
                        box.padding = 0.075,
                        fontface = "italic"
        ) +
        coord_cartesian(expand = T, clip = "off") +
        ylim(-1, 3)
        
      
    } else if (p2.labels == TRUE & length(geneSetID) == 1) {
      gsdata_labels <- gsdata %>% filter(ymax > 0)
      max_y <- max(gsdata_labels$ymax)
      
      description_labels <- unique(gsdata_labels$Description) 
      description_cut_point <- ceiling(length(description_labels) / 2)
      description_top <- description_labels[(description_cut_point+1):length(description_labels)]
      
      gsdata_labels_top <- gsdata_labels %>% filter(Description %in% description_top)
      max_y <- max(gsdata_labels$ymax)
      p2 <- p2 + 
        geom_text_repel(data = gsdata_labels_top,
                        aes(y = max_y, label = gene),
                        force_pull   = 0, # do not pull toward data points
                        nudge_y      = max_y * 0.05,
                        direction    = "x",
                        ylim = c(1.4, 2),
                        angle        = 90,
                        hjust        = 0,
                        max.time = 3,
                        max.overlaps = Inf,
                        min.segment.length = 0,
                        segment.size = 0.2,
                        segment.square = TRUE,
                        segment.inflect = FALSE,
                        segment.curvature = 1,
                        segment.angle = -1,
                        box.padding = 0.075,
                        fontface = "italic"
        )
    } else {
      p2 <- p2 + theme(panel.background = element_rect(fill = p2.fill))
    }
    if (length(geneSetID) == 1 & p2.add.intervals.colors) {
        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) 
            inv <- inv + 1
        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
            xmax = xmax, col = col[unique(inv)])
        p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
            ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
            alpha = 0.9, inherit.aes = FALSE)
    }
    
    ##
    ## Rank plot
    ##
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + 
      geom_segment(data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + 
      ylab(p.pos.ylab) + 
      xlab("Rank in Ordered Dataset") + 
      theme(
        plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm"),
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill ="transparent",color = "black", linewidth = axis.linewidth)
        )
    
    ##
    ## Combine plots
    ##
    if (!is.null(title) && !is.na(title) && title != "") 
        p.res <- p.res + ggtitle(title)
    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }
        else if (p2.segment.color.var == 'Description') {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
    if (p2.segment.color.var %in% c('geneList','geneListCapped') ) {
      p2 <- p2 + scale_color_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c", midpoint = 0)
    }
    if (pvalue_table) {
        pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
        rownames(pd) <- pd$Description
        pd <- pd[, -1]
        for (i in seq_len(ncol(pd))) {
            pd[, i] <- format(pd[, i], digits = 4)
        }
        tp <- tableGrob2(pd, p.res)
        p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
            xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                0.95), ymin = quantile(p.res$data$runningScore, 
                0.75), ymax = quantile(p.res$data$runningScore, 
                0.9))
    }
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + 
      theme(
        axis.ticks.x = element_line(),
        axis.text.x = element_text()
        )
    if (length(subplots) == 1) 
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
            r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    if (length(rel_heights) > length(subplots)) 
        rel_heights <- rel_heights[subplots]
    aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}
```

### Supplementary Figure 7L

Gene set enrichment analysis of Amhr2 KO versus control RNA-seq datasets (as reported in Fig. 7F) interrogating the enrichment of a mammary cancer proliferative signature, comparing Py2T cells treated as indicated; NES, normalized enrichment score.


``` r
i <- "InVivo_AMHR2_KO–over–InVivo_Control_sgRNA"
use_gsets <- c('STAT1 target genes and interactors', 'SMAD target genes and interactors')

gse_manual_list[[i]]@result %>% 
  filter(ID %in% use_gsets) %>% 
  mutate(
    NES = round(NES, 2),
    pvalue = format.pval(pvalue, digits = 2),
    p.adjust = format.pval(p.adjust, digits = 2)
  ) %>% 
  dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
  kbl(caption = 'GSEA results',
      row.names = F) %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
              full_width = F) %>%
  column_spec(1, bold = TRUE) %>% 
  print
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>GSEA results</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:right;"> setSize </th>
   <th style="text-align:right;"> NES </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> p.adjust </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> STAT1 target genes and interactors </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 2.23 </td>
   <td style="text-align:left;"> 3.5e-06 </td>
   <td style="text-align:left;"> 7.0e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> SMAD target genes and interactors </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> -2.13 </td>
   <td style="text-align:left;"> 5.5e-06 </td>
   <td style="text-align:left;"> 7.4e-06 </td>
  </tr>
</tbody>
</table>

``` r
cat("\n\n")
```

``` r
# n_gsets <- seq_along(gse_manual_list[[i]]$Description)
n_gsets <- which(gse_manual_list[[i]]$Description %in% use_gsets)
for(j in n_gsets) {
  es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
  cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
  res <- gseaplot2_custom(
    gse_manual_list[[i]], 
    geneSetID = j, 
    title = i,
    base_size = 8,
    rel_heights = c(1, 1.5, 0),
    color =es_color,
    p2.fill = 'grey90',
    p2.linerange.linewidth = 0.80,
    legend.position = 'top',
    p.pos.ylab = 'Log2 Fold Change',
    p2.labels = TRUE,
    reverse.rank = FALSE,
    p2.segment.color.var = 'geneList',
    geneList.limit = 5
  )
  print(res)
  cat("\n\n")
}
```

#### STAT1 target genes and interactors 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7L-1.png" style="display: block; margin: auto;" />

#### SMAD target genes and interactors 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7L-2.png" style="display: block; margin: auto;" />


### Supplementary Figure 7M (Galunisertib positive)

Gene set enrichment analysis of Amhr2 KO versus control RNA-seq datasets (as reported in Fig. 7F) interrogating the enrichment of a mammary cancer proliferative signature, comparing Py2T cells treated in vitro with Galunisertib as indicated. 


``` r
i <- "AMHR2_KO_noTGFbeta--over--AMHR2_Control_culture"
use_gsets <- c('Mammary epithelium proliferation')

gse_manual_list[[i]]@result %>% 
  filter(ID %in% use_gsets) %>% 
  mutate(
    NES = round(NES, 2),
    pvalue = format.pval(pvalue, digits = 2),
    p.adjust = format.pval(p.adjust, digits = 2)
  ) %>% 
  dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
  kbl(caption = 'GSEA results',
      row.names = F) %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
              full_width = F) %>%
  column_spec(1, bold = TRUE) %>% 
  print
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>GSEA results</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:right;"> setSize </th>
   <th style="text-align:right;"> NES </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> p.adjust </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> Mammary epithelium proliferation </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> -0.65 </td>
   <td style="text-align:left;"> 0.9 </td>
   <td style="text-align:left;"> 0.99 </td>
  </tr>
</tbody>
</table>

``` r
cat("\n\n")
```

``` r
n_gsets <- which(gse_manual_list[[i]]$Description %in% use_gsets)
for(j in n_gsets) {
  es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
  cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
  res <- gseaplot2_custom(
    gse_manual_list[[i]], 
    geneSetID = j, 
    title = i,
    base_size = 8,
    rel_heights = c(1, 1.5, 0),
    color =es_color,
    p2.fill = 'grey90',
    p2.linerange.linewidth = 0.80,
    legend.position = 'top',
    p.pos.ylab = 'Log2 Fold Change',
    p2.labels = TRUE,
    reverse.rank = FALSE,
    p2.segment.color.var = 'geneList',
    geneList.limit = 5
  )
  print(res)
  cat("\n\n")
}
```

#### Mammary epithelium proliferation 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7Mtop-1.png" style="display: block; margin: auto;" />




### Supplementary Figure 7P

Gene set enrichment analysis of Amhr2 KO versus control RNA-seq datasets (as reported in Fig. 7F) interrogating the enrichment of a mammary cancer proliferative signature or the enrichment of an EMT signature including epithelial and mesenchymal genes, showing data from orthotopic Py2T tumors in FVB/N mice.


``` r
i <- "InVivo_AMHR2_KO + AMHR2_KO_culture--over--InVivo_Control_sgRNA + AMHR2_Control_culture"
use_gsets <- c('Mammary epithelium proliferation', 'EMT regulators and targets')

gse_manual_list[[i]]@result %>% 
  filter(ID %in% use_gsets) %>% 
  mutate(
    NES = round(NES, 2),
    pvalue = format.pval(pvalue, digits = 2),
    p.adjust = format.pval(p.adjust, digits = 2)
  ) %>% 
  dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
  kbl(caption = 'GSEA results',
      row.names = F) %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
              full_width = F) %>%
  column_spec(1, bold = TRUE) %>% 
  print
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>GSEA results</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:right;"> setSize </th>
   <th style="text-align:right;"> NES </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> p.adjust </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> EMT regulators and targets </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> -2.38 </td>
   <td style="text-align:left;"> 2.1e-08 </td>
   <td style="text-align:left;"> 8.6e-08 </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;"> Mammary epithelium proliferation </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> -1.35 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.12 </td>
  </tr>
</tbody>
</table>

``` r
cat("\n\n")
```

``` r
n_gsets <- which(gse_manual_list[[i]]$Description %in% use_gsets)
for(j in n_gsets) {
  es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
  cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
  res <- gseaplot2_custom(
    gse_manual_list[[i]], 
    geneSetID = j, 
    title = i,
    base_size = 8,
    rel_heights = c(1, 1.5, 0),
    color =es_color,
    p2.fill = 'grey90',
    p2.linerange.linewidth = 0.80,
    legend.position = 'top',
    p.pos.ylab = 'Log2 Fold Change',
    p2.labels = TRUE,
    reverse.rank = FALSE,
    p2.segment.color.var = 'geneList',
    geneList.limit = 5
  )
  print(res)
  cat("\n\n")
}
```

#### EMT regulators and targets 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7p-1.png" style="display: block; margin: auto;" />

#### Mammary epithelium proliferation 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7p-2.png" style="display: block; margin: auto;" />

### Supplementary Figure 7U left (Galunisertib negative)

Gene set enrichment analysis of Amhr2 KO versus control RNA-seq datasets (as reported in Fig. 7F) interrogating the enrichment of an IFN-α transcriptional target signature, comparing Py2T cells treated in vitro with Galunisertib as indicated.


``` r
i <- "InVivo_AMHR2_KO + AMHR2_KO_culture--over--InVivo_Control_sgRNA + AMHR2_Control_culture"
use_gsets <- c('STAT1 target genes and interactors')

gse_manual_list[[i]]@result %>% 
  filter(ID %in% use_gsets) %>% 
  mutate(
    NES = round(NES, 2),
    pvalue = format.pval(pvalue, digits = 2),
    p.adjust = format.pval(p.adjust, digits = 2)
  ) %>% 
  dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
  kbl(caption = 'GSEA results',
      row.names = F) %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
              full_width = F) %>%
  column_spec(1, bold = TRUE) %>% 
  print
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>GSEA results</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:right;"> setSize </th>
   <th style="text-align:right;"> NES </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> p.adjust </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> STAT1 target genes and interactors </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 1.97 </td>
   <td style="text-align:left;"> 0.00014 </td>
   <td style="text-align:left;"> 0.00027 </td>
  </tr>
</tbody>
</table>

``` r
cat("\n\n")
```

``` r
n_gsets <- which(gse_manual_list[[i]]$Description %in% use_gsets)
for(j in n_gsets) {
  es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
  cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
  res <- gseaplot2_custom(
    gse_manual_list[[i]], 
    geneSetID = j, 
    title = i,
    base_size = 8,
    rel_heights = c(1, 1.5, 0),
    color =es_color,
    p2.fill = 'grey90',
    p2.linerange.linewidth = 0.80,
    legend.position = 'top',
    p.pos.ylab = 'Log2 Fold Change',
    p2.labels = TRUE,
    reverse.rank = FALSE,
    p2.segment.color.var = 'geneList',
    geneList.limit = 5
  )
  print(res)
  cat("\n\n")
}
```

#### STAT1 target genes and interactors 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7Uleft-1.png" style="display: block; margin: auto;" />


### Supplementary Figure 7U right (Galunisertib positive)

Gene set enrichment analysis of Amhr2 KO versus control RNA-seq datasets (as reported in Fig. 7F) interrogating the enrichment of an IFN-α transcriptional target signature, comparing Py2T cells treated in vitro with Galunisertib as indicated.


``` r
i <- "AMHR2_KO_noTGFbeta--over--AMHR2_Control_culture"
use_gsets <- c('STAT1 target genes and interactors')

gse_manual_list[[i]]@result %>% 
  filter(ID %in% use_gsets) %>% 
  mutate(
    NES = round(NES, 2),
    pvalue = format.pval(pvalue, digits = 2),
    p.adjust = format.pval(p.adjust, digits = 2)
  ) %>% 
  dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
  kbl(caption = 'GSEA results',
      row.names = F) %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
              full_width = F) %>%
  column_spec(1, bold = TRUE) %>% 
  print
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>GSEA results</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:right;"> setSize </th>
   <th style="text-align:right;"> NES </th>
   <th style="text-align:left;"> pvalue </th>
   <th style="text-align:left;"> p.adjust </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;"> STAT1 target genes and interactors </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 1.4 </td>
   <td style="text-align:left;"> 0.079 </td>
   <td style="text-align:left;"> 0.32 </td>
  </tr>
</tbody>
</table>

``` r
cat("\n\n")
```

``` r
n_gsets <- which(gse_manual_list[[i]]$Description %in% use_gsets)
for(j in n_gsets) {
  es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
  cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
  res <- gseaplot2_custom(
    gse_manual_list[[i]], 
    geneSetID = j, 
    title = i,
    base_size = 8,
    rel_heights = c(1, 1.5, 0),
    color =es_color,
    p2.fill = 'grey90',
    p2.linerange.linewidth = 0.80,
    legend.position = 'top',
    p.pos.ylab = 'Log2 Fold Change',
    p2.labels = TRUE,
    reverse.rank = FALSE,
    p2.segment.color.var = 'geneList',
    geneList.limit = 5
  )
  print(res)
  cat("\n\n")
}
```

#### STAT1 target genes and interactors 

<img src="figure/rnaseq-amhr2_ko-deg.Rmd/gseaplot2-custom-labels-individual-plot-sf7Uright-1.png" style="display: block; margin: auto;" />

## All GSEA plots

``` r
i <- names(gse_manual_list)[1]
for(i in names(gse_manual_list)) {
  cat("### Comparison", i, "\n\n")
  
  gse_manual_list[[i]]@result %>% 
    mutate(
      NES = round(NES, 2),
      pvalue = format.pval(pvalue, digits = 2),
      p.adjust = format.pval(p.adjust, digits = 2)
    ) %>% 
    dplyr::select(Description, setSize, NES, pvalue, p.adjust) %>% 
    kbl(caption = 'GSEA results',
        row.names = F) %>%
    kable_paper(bootstrap_options = c("striped", "hover", "condensed"), 
                full_width = F) %>%
    column_spec(1, bold = TRUE) %>% 
    print
  
  cat("\n\n")
  
  n_gsets <- seq_along(gse_manual_list[[i]]$Description)
  for(j in n_gsets) {
    es_color <- ifelse(gse_manual_list[[i]]$NES[j] > 0, "#e41a1c", "#377eb8")
    cat("####", gse_manual_list[[i]]$Description[j], "\n\n")
    res <- gseaplot2_custom(
      gse_manual_list[[i]], 
      geneSetID = j, 
      title = i,
      base_size = 8,
      rel_heights = c(1, 1.5, 0),
      color =es_color,
      p2.fill = 'grey90',
      p2.linerange.linewidth = 0.80,
      legend.position = 'top',
      p.pos.ylab = 'Log2 Fold Change',
      p2.labels = TRUE,
      reverse.rank = FALSE,
      p2.segment.color.var = 'geneList',
      geneList.limit = 5
    )
    print(res)
    cat("\n\n")
  }
}
```

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
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] circlize_0.4.16             ComplexHeatmap_2.22.0      
 [3] org.Mm.eg.db_3.20.0         AnnotationDbi_1.68.0       
 [5] PCAtools_2.18.0             pheatmap_1.0.13            
 [7] scran_1.34.0                scuttle_1.16.0             
 [9] SingleCellExperiment_1.28.1 DESeq2_1.46.0              
[11] edgeR_4.4.2                 limma_3.62.2               
[13] SummarizedExperiment_1.36.0 Biobase_2.66.0             
[15] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[17] IRanges_2.40.1              S4Vectors_0.44.0           
[19] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
[21] matrixStats_1.5.0           openxlsx_4.2.8             
[23] magrittr_2.0.3              ggrepel_0.9.6              
[25] ggbeeswarm_0.7.2            RColorBrewer_1.1-3         
[27] colorblindr_0.1.0           colorspace_2.1-1           
[29] cowplot_1.2.0               kableExtra_1.4.0           
[31] knitr_1.50                  DT_0.34.0                  
[33] foreach_1.5.2               showtext_0.9-7             
[35] showtextdb_3.0              sysfonts_0.8.9             
[37] lubridate_1.9.4             forcats_1.0.0              
[39] stringr_1.5.1               dplyr_1.1.4                
[41] purrr_1.1.0                 readr_2.1.5                
[43] tidyr_1.3.1                 tibble_3.3.0               
[45] ggplot2_3.5.2               tidyverse_2.0.0            
[47] workflowr_1.7.1            

loaded via a namespace (and not attached):
  [1] splines_4.4.3             later_1.4.4              
  [3] ggplotify_0.1.2           R.oo_1.27.1              
  [5] lifecycle_1.0.4           doParallel_1.0.17        
  [7] rprojroot_2.1.0           processx_3.8.6           
  [9] lattice_0.22-7            crosstalk_1.2.2          
 [11] sass_0.4.10               rmarkdown_2.29           
 [13] jquerylib_0.1.4           yaml_2.3.10              
 [15] ggtangle_0.0.7            metapod_1.14.0           
 [17] httpuv_1.6.16             zip_2.3.3                
 [19] DBI_1.2.3                 abind_1.4-8              
 [21] zlibbioc_1.52.0           R.utils_2.13.0           
 [23] yulab.utils_0.2.0         git2r_0.36.2             
 [25] GenomeInfoDbData_1.2.13   enrichplot_1.26.6        
 [27] irlba_2.3.5.1             tidytree_0.4.6           
 [29] dqrng_0.4.1               svglite_2.2.1            
 [31] DelayedMatrixStats_1.28.1 codetools_0.2-20         
 [33] DelayedArray_0.32.0       DOSE_4.0.1               
 [35] xml2_1.3.8                tidyselect_1.2.1         
 [37] shape_1.4.6.1             aplot_0.2.8              
 [39] UCSC.utils_1.2.0          farver_2.1.2             
 [41] ScaledMatrix_1.14.0       jsonlite_2.0.0           
 [43] GetoptLong_1.0.5          BiocNeighbors_2.0.1      
 [45] iterators_1.0.14          systemfonts_1.2.3        
 [47] tools_4.4.3               treeio_1.30.0            
 [49] Rcpp_1.1.0                glue_1.8.0               
 [51] SparseArray_1.6.2         xfun_0.53                
 [53] here_1.0.1                qvalue_2.38.0            
 [55] withr_3.0.2               fastmap_1.2.0            
 [57] bluster_1.16.0            callr_3.7.6              
 [59] digest_0.6.37             rsvd_1.0.5               
 [61] gridGraphics_0.5-1        timechange_0.3.0         
 [63] R6_2.6.1                  textshaping_1.0.1        
 [65] GO.db_3.20.0              dichromat_2.0-0.1        
 [67] RSQLite_2.4.3             R.methodsS3_1.8.2        
 [69] generics_0.1.4            data.table_1.17.8        
 [71] httr_1.4.7                htmlwidgets_1.6.4        
 [73] S4Arrays_1.6.0            whisker_0.4.1            
 [75] pkgconfig_2.0.3           gtable_0.3.6             
 [77] blob_1.2.4                XVector_0.46.0           
 [79] htmltools_0.5.8.1         fgsea_1.32.4             
 [81] clue_0.3-66               scales_1.4.0             
 [83] png_0.1-8                 ggfun_0.2.0              
 [85] rstudioapi_0.17.1         tzdb_0.5.0               
 [87] reshape2_1.4.4            rjson_0.2.23             
 [89] nlme_3.1-168              cachem_1.1.0             
 [91] GlobalOptions_0.1.2       parallel_4.4.3           
 [93] vipor_0.4.7               pillar_1.11.0            
 [95] vctrs_0.6.5               promises_1.3.3           
 [97] BiocSingular_1.22.0       beachmat_2.22.0          
 [99] cluster_2.1.8.1           beeswarm_0.4.0           
[101] evaluate_1.0.5            cli_3.6.5                
[103] locfit_1.5-9.12           compiler_4.4.3           
[105] rlang_1.1.6               crayon_1.5.3             
[107] labeling_0.4.3            ps_1.9.1                 
[109] getPass_0.2-4             plyr_1.8.9               
[111] fs_1.6.6                  stringi_1.8.7            
[113] viridisLite_0.4.2         BiocParallel_1.40.2      
[115] Biostrings_2.74.1         lazyeval_0.2.2           
[117] GOSemSim_2.32.0           Matrix_1.7-3             
[119] patchwork_1.3.1           hms_1.1.3                
[121] sparseMatrixStats_1.18.0  bit64_4.6.0-1            
[123] KEGGREST_1.46.0           statmod_1.5.0            
[125] igraph_2.1.4              memoise_2.0.1            
[127] bslib_0.9.0               ggtree_3.14.0            
[129] fastmatch_1.1-6           bit_4.6.0                
[131] ape_5.8-1                
</div>
