---
title: "StealTHY CRISPR screen"
subtitle: "Mouse libray (n = 2215)"
author: "Francesc Castro-Giner"
date: 'October 31, 2025'
output: workflowr::wflow_html
editor_options:
  chunk_output_type: console
params:
  data_dir: ./data/crispr/mm_2215_sgRNA
  se_rnaseq_bulk: ./data/rnaseq/bulk/se_gene.rds
  se_rnaseq_facs: ./data/rnaseq/facs/se_qc.rds
  output_dir: ./output/crispr/mm_2215_sgRNA
  mageck_dir: mageck_rra_expr
  ne_mageck_dir: mageck_rra
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
	Untracked:  analysis/crispr-mm_2215_sgRNA-clonality.md
	Untracked:  analysis/index.md
	Untracked:  analysis/rnaseq-tumor-bulk.md
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
to the R Markdown (<code>analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd</code>) and HTML (<code>docs/crispr-mm_2215_sgRNA-StealTHY.html</code>)
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
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/455fa9121fcca51363cf639f7ee170fb108e34c1/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">455fa91</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Prepare files for publication</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/crispr-mm_2215_sgRNA-StealTHY.html" target="_blank">32aca02</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/e391b94d2101a3ebbb4ef002cf4887c5d14738d0/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">e391b94</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>update crispr plots</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/00ae28bcb2715ec6f4030d74d92f51359fa43ce4/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">00ae28b</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-05-05</td>
<td>updated r1oa</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/1c45f78451a8507c8de37da9a7d78c2e7cd92a3a/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">1c45f78</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-04-30</td>
<td>rm temporary files</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/5fab19d765a869233e820f2aed94595a89d5fe9b/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">5fab19d</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-28</td>
<td>clean code</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/1920712c9ede8023880e174fc87f3cd238c3d612/analysis/crispr-mm_2215_sgRNA-StealTHY.Rmd" target="_blank">1920712</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-24</td>
<td>added F5</td>
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

library(MAGeCKFlute)
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
  rownames

se <- se[,use_cols]
```

Data wrangling

``` r
chunk_colData <- colData(se) %>% data.frame %>%
  mutate(
    vector_short = ifelse(vector_short == 'Puro_EF1alpha_SpCas9', 
                          'Cas9+Puro', vector_short),
    vector_short = factor(vector_short, levels = c('Thy1', 'Puro', 'Cas9+Puro')),
    mouse_model_short = case_match(mouse_model_global,
                                   'immunodeficient' ~ 'NSG',
                                   'immunocompetent' ~ 'Syngeneic'
                                   ),
    mouse_model_short = factor(mouse_model_short, levels = c('NSG', 'Syngeneic')),
    vector_mhost = paste0(vector_short, ', ', mouse_model_short),
    model_vector_mhost = paste0(donor,', ', vector_short, ', ', mouse_model_short),
    
    sample_type = case_match(sample_type, 
                             'whole_tumor' ~ 'Primary tumor', 
                             'whole_lung' ~ 'Lung mets'),
    sample_type = factor(sample_type, levels = c('Primary tumor', 'Lung mets'))
  )

colData(se) <- chunk_colData %>% DataFrame
```


## Configure analyses
### sgRNA distribution
Define comparisons for modality and Kolmogorov-Smirnov test

``` r
x <- colData(se) %>% data.frame

group_list_modality <- list(
  `StealTHY syngeneic` = x %>% 
    filter(vector_mmodel %in% c('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Primary tumor') %>%
    pull(sample_alias),
  
  `Cas9+Puro syngeneic` = x %>% 
    filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Primary tumor') %>% 
    pull(sample_alias),
  
  `StealTHY NSG` = x %>% 
    filter(vector_mmodel %in% c('4T1-KO-Thy1-NSG', 'MVT1-KO-Thy1-NSG', 'Py2T-KO-Thy1-NSG') & sample_type == 'Primary tumor') %>%
    pull(sample_alias),
  
  `Cas9+Puro NSG` = x %>% 
    filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-NSG' & sample_type == 'Primary tumor') %>% 
    pull(sample_alias),
  
  `StealTHY syngeneic 4T1` = x %>% 
    filter(vector_mmodel == '4T1-KO-Thy1-BALB' & sample_type == 'Primary tumor') %>%
    pull(sample_alias)
  
)


comp_list_modality <- list(
  `StealTHY syngeneic vs. StealTHY NSG` = list(
    `StealTHY syngeneic` = group_list_modality$`StealTHY syngeneic`,
    `StealTHY NSG` = group_list_modality$`StealTHY NSG`
  ),
  `StealTHY syngeneic vs. Cas9+Puro` = list(
    `StealTHY syngeneic` = group_list_modality$`StealTHY syngeneic 4T1`,
    `Cas9+Puro syngeneic` = group_list_modality$`Cas9+Puro syngeneic`
  )
)

group_df_modality <- foreach(i = names(group_list_modality), .combine = rbind) %do% {
  data.frame(
    sample_alias = group_list_modality[[i]],
    modality_group = i
  )
}
```

### MAGeCK analysis

#### List of comparisons
List of comparisons for differential abundance using MAGeCK

``` r
x <- colData(se) %>% data.frame

mageck_comp_list <- list(
  `StealTHY KO Syngeneic` = list(
    `Primary tumor` = x %>% filter(vector_mmodel %in% list('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel %in% list('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  `StealTHY KO NSG` = list(
    `Primary tumor` = x %>% filter(global_vector_mmodel == 'whole_tumor-KO-Thy1-immunodeficient') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(global_vector_mmodel == 'whole_lung-KO-Thy1-immunodeficient') %>% pull(sample_alias)
  ),
  `Cas9+Puro Syngeneic` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Lung mets') %>% pull(sample_alias)
  ),
  
  `StealTHY KO over dCas9 In vitro` = list(
    `StealTHY` = x %>% filter(global_vector_mmodel_invitro_cas9transfection == 'bulk-KO-Thy1-Invitro-pre-Cas9') %>% pull(sample_alias),
    `dCas9` = x %>% filter(global_vector_mmodel_invitro_cas9transfection_grouped == 'KO-Thy1-Invitro-post-Cas9') %>% pull(sample_alias)
    ),
  
  `StealTHY KO over dCas9 In vivo` = list(
    `StealTHY` = x %>% filter(vector_mmodel %in% list('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `dCas9 mets` = x %>% filter(vector_mmodel %in% list('4T1-noKO-Thy1-BALB', 'MVT1-noKO-Thy1-FVB', 'Py2T-noKO-Thy1-FVB') & sample_type == 'Primary tumor') %>% pull(sample_alias)
    )
)


mageck_comp_list_names <- names(mageck_comp_list) 
names(mageck_comp_list_names) <- names(mageck_comp_list) %>% gsub(" ", "_",.) %>% tolower
```


``` r
data_comp <- foreach(i = names(mageck_comp_list), .combine = rbind) %do% {
  c(
    Comparison = i,
    Case = paste(mageck_comp_list[[i]][[1]], collapse = " "),
    Control = paste(mageck_comp_list[[i]][[2]], collapse = " ")
  )
}

data_comp %>% 
  datatable(., 
            rownames = FALSE, 
            filter = 'top', 
            caption = 'List of comparisons for MAGeCK analysis',
            extensions = 'Buttons', 
            options = list(
              dom = 'Blfrtip',
              buttons = c('csv', 'excel'),
              columnDefs = list(
                list(width = 120, targets = 1:2)
              )
            )
            )
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-0e41329c1963262c826a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-0e41329c1963262c826a">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>List of comparisons for MAGeCK analysis<\/caption>","data":[["StealTHY KO Syngeneic","StealTHY KO NSG","Cas9+Puro Syngeneic","StealTHY KO over dCas9 In vitro","StealTHY KO over dCas9 In vivo"],["CRISPRKO_Clone_1 CRISPRKO_Clone_10 CRISPRKO_Clone_11 CRISPRKO_Clone_12 CRISPRKO_Clone_13 CRISPRKO_Clone_14 CRISPRKO_Clone_2 CRISPRKO_Clone_20 CRISPRKO_Clone_21 CRISPRKO_Clone_22 CRISPRKO_Clone_23 CRISPRKO_Clone_24 CRISPRKO_Clone_3 CRISPRKO_Clone_4 CRISPRKO","CRISPRKO_Clone_45 CRISPRKO_Clone_46 CRISPRKO_Clone_47 CRISPRKO_Clone_48 CRISPRKO_Clone_49 CRISPRKO_Clone_5 CRISPRKO_Clone_50 CRISPRKO_Clone_51 CRISPRKO_Clone_52 CRISPRKO_Clone_53 CRISPRKO_Clone_54 CRISPRKO_Clone_6 CRISPRKO_Clone_7 CRISPRKO_Clone_8 CRISPRKO_Clone_9","PlatePositionA5 PlatePositionA6 PlatePositionC5 PlatePositionC6 PlatePositionE5 PlatePositionE6 PlatePositionG4 PlatePositionG5 PlatePositionG6","CRISPRKO_Clone_61 CRISPRKO_Clone_64 CRISPRKO_Clone_67","CRISPRKO_Clone_1 CRISPRKO_Clone_10 CRISPRKO_Clone_11 CRISPRKO_Clone_12 CRISPRKO_Clone_13 CRISPRKO_Clone_14 CRISPRKO_Clone_2 CRISPRKO_Clone_20 CRISPRKO_Clone_21 CRISPRKO_Clone_22 CRISPRKO_Clone_23 CRISPRKO_Clone_24 CRISPRKO_Clone_3 CRISPRKO_Clone_4 CRISPRKO"],["CRISPRKO_Clone_15 CRISPRKO_Clone_16 CRISPRKO_Clone_17 CRISPRKO_Clone_18 CRISPRKO_Clone_19 CRISPRKO_Clone_30 CRISPRKO_Clone_31 CRISPRKO_Clone_32 CRISPRKO_Clone_33 CRISPRKO_Clone_34 CRISPRKO_Clone_35 CRISPRKO_Clone_36 CRISPRKO_Clone_37 CRISPRKO_Clone_38 CRISPRKO_Clone_39","CRISPRKO_Clone_25 CRISPRKO_Clone_26 CRISPRKO_Clone_27 CRISPRKO_Clone_28 CRISPRKO_Clone_29 CRISPRKO_Clone_40 CRISPRKO_Clone_41 CRISPRKO_Clone_42 CRISPRKO_Clone_43 CRISPRKO_Clone_44 CRISPRKO_Clone_55 CRISPRKO_Clone_56 CRISPRKO_Clone_57 CRISPRKO_Clone_58 CRISPRKO_Clone_59","PlatePositionB5 PlatePositionD5 PlatePositionD6 PlatePositionF3 PlatePositionF5 PlatePositionF6 PlatePositionH4 PlatePositionH5 PlatePositionH6","CRISPRKO_Clone_60 CRISPRKO_Clone_65 CRISPRKO_Clone_68 PlatePositionA11 PlatePositionB11 PlatePositionC11 PlatePositionD11 PlatePositionG10 PlatePositionH10","4T_TumoT_B1 4T_TumoT_B2 4T_TumoT_B3 4T_TumoT_B4 4T_TumoT_B5 MVT_TumoT_F1 MVT_TumoT_F2 MVT_TumoT_F3 MVT_TumoT_F4 MVT_TumoT_F5 PY_TumoT_F1 PY_TumoT_F2 PY_TumoT_F3 PY_TumoT_F4 PY_TumoT_F5"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Comparison<\/th>\n      <th>Case<\/th>\n      <th>Control<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"width":120,"targets":[1,2]},{"name":"Comparison","targets":0},{"name":"Case","targets":1},{"name":"Control","targets":2}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>
```

#### Genes to remove

``` r
min_counts <- 1
min_prop <- 0.30
```

Remove genes not expressed in the RNA-seq experiment: we keep genes that have counts above 1 counts in 30 % of the samples from at least one group. The group is defined by the cell model. We consider only non-NSG samples for this analysis. additionally we also remove *Gt(ROSA)26Sor* and *Cd4* (they cannot be expressed by cancer cells).


``` r
# Read RNA-seq data
se_rna <- readRDS(params$se_rnaseq_bulk)
assay(se_rna, 'cpm') <- edgeR::cpm(se_rna)

# CRISPR library genes
crispr_genes <- unique(rowData(se)$Gene) %>% grep('Non_Target', ., invert = T, value = T)

# RNA-seq : Select genes and non-NSG samples
use_cols <- colData(se_rna)$mouse_model != 'NSG'
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_cols]

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[rowSums(counts_stats >= min_prop) > 0]

# Change to gene names
genes_to_keep <- rowData(use_se[rows_to_keep,])$gene_name %>% 
  unique

# Remove Rosa26 and Cd4 : cannot be expressed by cancer cells
genes_to_keep <- genes_to_keep[!genes_to_keep %in% c('Gt(ROSA)26Sor', 'Cd4')]

# Define genes to remove
genes_to_rm <- setdiff(crispr_genes, genes_to_keep)

# Show stats for removed genes
counts_stats$gene_name <- rowData(use_se[rownames(counts_stats),])$gene_name
counts_stats %<>% dplyr::select(gene_name, everything()) 

counts_stats_filterByExpr <- counts_stats

# Filter Crispr Data
keep_gene_guides <- rownames(se)[rowData(se)$Gene %in% genes_to_keep]
keep_ctrl_guides <- rownames(se)[grepl('Non_Target', rowData(se)$Gene)]
```

The table below shows the proportion of samples with counts >= 1in each group for genes removed (n = 66`).

``` r
counts_stats_filterByExpr %>% 
  filter(gene_name %in% genes_to_rm) %>% 
  arrange(desc(all)) %>% 
  datatable(., 
            rownames = FALSE, 
            filter = 'top', 
            caption = 'Statistics',
            extensions = 'Buttons', 
            options = list(
              dom = 'Blfrtip',
              buttons = c('csv', 'excel')
              )) %>% 
  formatRound(columns = 2:5)
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-cf74d9caecbce7389320" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-cf74d9caecbce7389320">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"1\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.916666666666667\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.416666666666667\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"0.777777777777778\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>Statistics<\/caption>","data":[["Cd4","Bmp5","Acvr1c","Il12a","Il4","Mpl","Fgf23","Gdf2","Cxcl15","Il31","Lep","Ccr1l1","Nrg3","Cd40lg","Il21","Angptl3","Il5ra","Fgf3","Il13","Fgf14","Tnfrsf17","Cd70","Fgf8","Cxcr1","Ifnk","Ifna12","Ifna2","Fgf6","Fgf4","Nodal","Il9","Il17a","Mstn","Erbb4","Il1f10","Ccl21b","Ccl27b","Il11ra2","Ifna15","Ifna14","Ifna9","Ifna7","Ifna11","Ifna6","Ifna5","Ifna4","Ifna1","Ccl26","Ifnl2","Ifnl3","Klk1b4","Inhbc","Igfbp1","Il3","Prl","Il25","Fgf12","Dkk1","Fgf16"],[1,0.25,0.25,0.25,0.25,0.08333333333333333,0.25,0.08333333333333333,0.1666666666666667,0.1666666666666667,0,0,0,0.1666666666666667,0.1666666666666667,0,0.08333333333333333,0,0.1666666666666667,0,0.1666666666666667,0,0.08333333333333333,0,0.08333333333333333,0,0,0,0,0,0.08333333333333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.9166666666666666,0.25,0.08333333333333333,0.1666666666666667,0.08333333333333333,0,0.08333333333333333,0.25,0,0,0.08333333333333333,0.1666666666666667,0,0.08333333333333333,0,0.1666666666666667,0.08333333333333333,0,0,0.08333333333333333,0,0.1666666666666667,0,0.08333333333333333,0,0.08333333333333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.4166666666666667,0.1666666666666667,0.08333333333333333,0,0.08333333333333333,0.25,0,0,0.08333333333333333,0.08333333333333333,0.1666666666666667,0.08333333333333333,0.25,0,0,0,0,0.1666666666666667,0,0.08333333333333333,0,0,0.08333333333333333,0,0,0,0.08333333333333333,0.08333333333333333,0.08333333333333333,0.08333333333333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.7777777777777778,0.2222222222222222,0.1388888888888889,0.1388888888888889,0.1388888888888889,0.1111111111111111,0.1111111111111111,0.1111111111111111,0.08333333333333333,0.08333333333333333,0.08333333333333333,0.08333333333333333,0.08333333333333333,0.08333333333333333,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.05555555555555555,0.02777777777777778,0.02777777777777778,0.02777777777777778,0.02777777777777778,0.02777777777777778,0.02777777777777778,0.02777777777777778,0.02777777777777778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>gene_name<\/th>\n      <th>4T1<\/th>\n      <th>MVT1<\/th>\n      <th>Py2T<\/th>\n      <th>all<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":2,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":3,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"targets":4,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 2, 3, \",\", \".\", null);\n  }"},{"className":"dt-right","targets":[1,2,3,4]},{"name":"gene_name","targets":0},{"name":"4T1","targets":1},{"name":"MVT1","targets":2},{"name":"Py2T","targets":3},{"name":"all","targets":4}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render"],"jsHooks":[]}</script>
```

#### Run MAGeCK

Generate MAGeCK test scripts. To generate the script files change `eval = FALSE` to `eval = TRUE` in the chunk options.

The scripts must be run in the terminal, inside the path ./output/crispr/mm_2215_sgRNA/mageck_rra_expr. The script will activate the mageck environment and run the test for each comparison. The results will be saved in the same directory.

``` r
chunk_se <- se[c(keep_gene_guides, keep_ctrl_guides),]
out_dir <- file.path(params$output_dir, params$mageck_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Generate count matrix
count_mat <- assay(chunk_se, 'counts') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('sgRNA') %>% 
  mutate(Gene = rowData(chunk_se)$Gene) %>% 
  dplyr::select(sgRNA, Gene, everything())
write_tsv(count_mat, file = file.path(out_dir, 'counts.txt'))

# Generate list of control sgRNA IDs
non_target_list <- rownames(chunk_se) %>% grep('Non_Target', ., value = T) %>% data.frame
write_tsv(non_target_list, file = file.path(out_dir, 'control_sgrna.txt'), col_names = F)

# Generate mageck commands
res_cmd <- data.frame('source activate mageckenv')
res_cmd <- foreach(i=names(mageck_comp_list), .combine = rbind) %do% {
  paste('mageck test -k counts.txt --control-sgrna control_sgrna.txt --norm-method control -t', paste(mageck_comp_list[[i]][[1]], collapse = ','),'-c',  paste(mageck_comp_list[[i]][[2]], collapse = ','), '-n', tolower(gsub(" ", "_", i)), '\n')
} %>% data.frame()

res_cmd <- rbind('conda activate mageckenv', res_cmd, 'conda deactivate')
write_tsv(res_cmd, file.path(out_dir, 'run_mageck_test.sh'), col_names = FALSE)
```

#### Load MAGeCK results

``` r
gene_summ_files <- list.files(path = file.path(params$output_dir, params$mageck_dir), pattern = 'gene_summary.txt', full.names = TRUE)
analysis_prefix <- basename(gene_summ_files) %>% gsub(".gene_summary.txt", "", .)
gene_summ <- foreach(i = gene_summ_files) %do% read.delim(i, check.names = FALSE)
names(gene_summ) <- analysis_prefix

sgrna_summ_files <- list.files(path = file.path(params$output_dir, params$mageck_dir), pattern = 'sgrna_summary.txt', full.names = TRUE)
analysis_prefix <- basename(sgrna_summ_files) %>% gsub(".sgrna_summary.txt", "", .)
sgrna_summ <- foreach(i = sgrna_summ_files) %do% read.delim(i, check.names = FALSE)
names(sgrna_summ) <- analysis_prefix
```


### MAGeCK analysis with non-expressed genes

#### List of comparisons
List of comparisons for differential abundance using MAGeCK

``` r
x <- colData(se) %>% data.frame

ne_mageck_comp_list <- list(
  
  `StealTHY` = list(
    `Primary tumor` = x %>% filter(vector_mmodel %in% list('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel %in% list('4T1-KO-Thy1-BALB', 'MVT1-KO-Thy1-FVB', 'Py2T-KO-Thy1-FVB') & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  `Puro` = list(
    `Primary tumor` = x %>% filter(global_vector_mmodel == 'whole_tumor-KO-Puro-immunocompetent' ) %>% pull(sample_alias),
    `Lung mets` = x %>% filter(global_vector_mmodel == 'whole_lung-KO-Puro-immunocompetent' ) %>% pull(sample_alias)
  ),
  
  `Cas9+Puro` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Lung mets') %>% pull(sample_alias)
  ),
  
  `StealTHY NSG` = list(
    `Primary tumor` = x %>% filter(global_vector_mmodel == 'whole_tumor-KO-Thy1-immunodeficient') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(global_vector_mmodel == 'whole_lung-KO-Thy1-immunodeficient') %>% pull(sample_alias)
  ),
  
  ## 4T1
  `StealTHY 4T1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Thy1-BALB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == '4T1-KO-Thy1-BALB' & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  `Puro 4T1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Puro-BALB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `In vitro` = x %>% filter(vector_mmodel == '4T1-KO-lentiguide-Puro-In vitro' & invitro_type == 'Puro-resistant') %>% pull(sample_alias)
  ),
  
  `Cas9+Puro 4T1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == '4T1-KO-Puro_EF1alpha_SpCas9-BALB' & sample_type == 'Lung mets') %>% pull(sample_alias)
  ),
  
  `StealTHY NSG 4T1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == '4T1-KO-Thy1-NSG' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == '4T1-KO-Thy1-NSG' & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  ## MVT1
  `StealTHY MVT1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'MVT1-KO-Thy1-FVB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == 'MVT1-KO-Thy1-FVB' & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  `Puro MVT1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'MVT1-KO-Puro-FVB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `In vitro` = x %>% filter(vector_mmodel == 'MVT1-KO-lentiguide-Puro-In vitro' & invitro_type == 'Puro-resistant') %>% pull(sample_alias)
  ),
  
  `StealTHY NSG MVT1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'MVT1-KO-Thy1-NSG' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == 'MVT1-KO-Thy1-NSG' & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  ## Py2T
  `StealTHY Py2T` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'Py2T-KO-Thy1-FVB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == 'Py2T-KO-Thy1-FVB' & sample_type == 'Lung mets') %>% pull(sample_alias)
    ),
  
  `Puro Py2T` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'Py2T-KO-Puro-FVB' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `In vitro` = x %>% filter(vector_mmodel == 'Py2T-KO-lentiguide-Puro-In vitro' & invitro_type == 'Puro-resistant') %>% pull(sample_alias)
  ),
  
  `StealTHY NSG Py2T` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'Py2T-KO-Thy1-NSG' & sample_type == 'Primary tumor') %>% pull(sample_alias),
    `Lung mets` = x %>% filter(vector_mmodel == 'Py2T-KO-Thy1-NSG' & sample_type == 'Lung mets') %>% pull(sample_alias)
    )
  
  
)


ne_mageck_comp_list_names <- names(ne_mageck_comp_list) 
names(ne_mageck_comp_list_names) <- names(ne_mageck_comp_list) %>% gsub(" ", "_",.) %>% tolower
```


``` r
data_comp <- foreach(i = names(ne_mageck_comp_list), .combine = rbind) %do% {
  c(
    Comparison = i,
    Case = paste(ne_mageck_comp_list[[i]][[1]], collapse = " "),
    Control = paste(ne_mageck_comp_list[[i]][[2]], collapse = " ")
  )
}

data_comp %>% 
  datatable(., 
            rownames = FALSE, 
            filter = 'top', 
            caption = 'List of comparisons for MAGeCK analysis',
            extensions = 'Buttons', 
            options = list(
              dom = 'Blfrtip',
              buttons = c('csv', 'excel'),
              columnDefs = list(
                list(width = 120, targets = 1:2)
              )
            )
            )
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-684d012f49e01d513276" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-684d012f49e01d513276">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>List of comparisons for MAGeCK analysis<\/caption>","data":[["StealTHY","Puro","Cas9+Puro","StealTHY NSG","StealTHY 4T1","Puro 4T1","Cas9+Puro 4T1","StealTHY NSG 4T1","StealTHY MVT1","Puro MVT1","StealTHY NSG MVT1","StealTHY Py2T","Puro Py2T","StealTHY NSG Py2T"],["CRISPRKO_Clone_1 CRISPRKO_Clone_10 CRISPRKO_Clone_11 CRISPRKO_Clone_12 CRISPRKO_Clone_13 CRISPRKO_Clone_14 CRISPRKO_Clone_2 CRISPRKO_Clone_20 CRISPRKO_Clone_21 CRISPRKO_Clone_22 CRISPRKO_Clone_23 CRISPRKO_Clone_24 CRISPRKO_Clone_3 CRISPRKO_Clone_4 CRISPRKO","PlatePositionA7 PlatePositionA8 PlatePositionA9 PlatePositionB7 PlatePositionB8 PlatePositionB9 PlatePositionC7 PlatePositionC8 PlatePositionC9 PlatePositionD7 PlatePositionD9","PlatePositionA5 PlatePositionA6 PlatePositionC5 PlatePositionC6 PlatePositionE5 PlatePositionE6 PlatePositionG4 PlatePositionG5 PlatePositionG6","CRISPRKO_Clone_45 CRISPRKO_Clone_46 CRISPRKO_Clone_47 CRISPRKO_Clone_48 CRISPRKO_Clone_49 CRISPRKO_Clone_5 CRISPRKO_Clone_50 CRISPRKO_Clone_51 CRISPRKO_Clone_52 CRISPRKO_Clone_53 CRISPRKO_Clone_54 CRISPRKO_Clone_6 CRISPRKO_Clone_7 CRISPRKO_Clone_8 CRISPRKO_Clone_9","CRISPRKO_Clone_10 CRISPRKO_Clone_11 CRISPRKO_Clone_12 CRISPRKO_Clone_13 CRISPRKO_Clone_14","PlatePositionA7 PlatePositionB7 PlatePositionC7 PlatePositionD7","PlatePositionA5 PlatePositionA6 PlatePositionC5 PlatePositionC6 PlatePositionE5 PlatePositionE6 PlatePositionG4 PlatePositionG5 PlatePositionG6","CRISPRKO_Clone_45 CRISPRKO_Clone_46 CRISPRKO_Clone_47 CRISPRKO_Clone_48 CRISPRKO_Clone_49","CRISPRKO_Clone_20 CRISPRKO_Clone_21 CRISPRKO_Clone_22 CRISPRKO_Clone_23 CRISPRKO_Clone_24","PlatePositionA8 PlatePositionB8 PlatePositionC8","CRISPRKO_Clone_5 CRISPRKO_Clone_6 CRISPRKO_Clone_7 CRISPRKO_Clone_8 CRISPRKO_Clone_9","CRISPRKO_Clone_1 CRISPRKO_Clone_2 CRISPRKO_Clone_3 CRISPRKO_Clone_4 CRISPRKO","PlatePositionA9 PlatePositionB9 PlatePositionC9 PlatePositionD9","CRISPRKO_Clone_50 CRISPRKO_Clone_51 CRISPRKO_Clone_52 CRISPRKO_Clone_53 CRISPRKO_Clone_54"],["CRISPRKO_Clone_15 CRISPRKO_Clone_16 CRISPRKO_Clone_17 CRISPRKO_Clone_18 CRISPRKO_Clone_19 CRISPRKO_Clone_30 CRISPRKO_Clone_31 CRISPRKO_Clone_32 CRISPRKO_Clone_33 CRISPRKO_Clone_34 CRISPRKO_Clone_35 CRISPRKO_Clone_36 CRISPRKO_Clone_37 CRISPRKO_Clone_38 CRISPRKO_Clone_39","PlatePositionE7 PlatePositionE8 PlatePositionE9 PlatePositionF7 PlatePositionF8 PlatePositionF9 PlatePositionG7 PlatePositionG8 PlatePositionG9 PlatePositionH7 PlatePositionH9","PlatePositionB5 PlatePositionD5 PlatePositionD6 PlatePositionF3 PlatePositionF5 PlatePositionF6 PlatePositionH4 PlatePositionH5 PlatePositionH6","CRISPRKO_Clone_25 CRISPRKO_Clone_26 CRISPRKO_Clone_27 CRISPRKO_Clone_28 CRISPRKO_Clone_29 CRISPRKO_Clone_40 CRISPRKO_Clone_41 CRISPRKO_Clone_42 CRISPRKO_Clone_43 CRISPRKO_Clone_44 CRISPRKO_Clone_55 CRISPRKO_Clone_56 CRISPRKO_Clone_57 CRISPRKO_Clone_58 CRISPRKO_Clone_59","CRISPRKO_Clone_15 CRISPRKO_Clone_16 CRISPRKO_Clone_17 CRISPRKO_Clone_18 CRISPRKO_Clone_19","CRISPRKO_Clone_71","PlatePositionB5 PlatePositionD5 PlatePositionD6 PlatePositionF3 PlatePositionF5 PlatePositionF6 PlatePositionH4 PlatePositionH5 PlatePositionH6","CRISPRKO_Clone_40 CRISPRKO_Clone_41 CRISPRKO_Clone_42 CRISPRKO_Clone_43 CRISPRKO_Clone_44","CRISPRKO_Clone_35 CRISPRKO_Clone_36 CRISPRKO_Clone_37 CRISPRKO_Clone_38 CRISPRKO_Clone_39","CRISPRKO_Clone_70","CRISPRKO_Clone_55 CRISPRKO_Clone_56 CRISPRKO_Clone_57 CRISPRKO_Clone_58 CRISPRKO_Clone_59","CRISPRKO_Clone_30 CRISPRKO_Clone_31 CRISPRKO_Clone_32 CRISPRKO_Clone_33 CRISPRKO_Clone_34","CRISPRKO_Clone_63","CRISPRKO_Clone_25 CRISPRKO_Clone_26 CRISPRKO_Clone_27 CRISPRKO_Clone_28 CRISPRKO_Clone_29"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Comparison<\/th>\n      <th>Case<\/th>\n      <th>Control<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"width":120,"targets":[1,2]},{"name":"Comparison","targets":0},{"name":"Case","targets":1},{"name":"Control","targets":2}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>
```


#### Run MAGeCK

Generate MAGeCK test scripts. To generate the script files change `eval = FALSE` to `eval = TRUE` in the chunk options.

The scripts must be run in the terminal, inside the path ./output/crispr/mm_2215_sgRNA/mageck_rra. The script will activate the mageck environment and run the test for each comparison. The results will be saved in the same directory.

``` r
chunk_se <- se
out_dir <- file.path(params$output_dir, params$ne_mageck_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Generate count matrix
count_mat <- assay(chunk_se, 'counts') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('sgRNA') %>% 
  mutate(Gene = rowData(chunk_se)$Gene) %>% 
  dplyr::select(sgRNA, Gene, everything())
write_tsv(count_mat, file = file.path(out_dir, 'counts.txt'))

# Generate list of control sgRNA IDs
non_target_list <- rownames(chunk_se) %>% grep('Non_Target', ., value = T) %>% data.frame
write_tsv(non_target_list, file = file.path(out_dir, 'control_sgrna.txt'), col_names = F)

# Generate mageck commands
res_cmd <- data.frame('source activate mageckenv')
res_cmd <- foreach(i=names(ne_mageck_comp_list), .combine = rbind) %do% {
  paste('mageck test -k counts.txt --control-sgrna control_sgrna.txt --norm-method control -t', paste(ne_mageck_comp_list[[i]][[1]], collapse = ','),'-c',  paste(ne_mageck_comp_list[[i]][[2]], collapse = ','), '-n', tolower(gsub(" ", "_", i)), '\n')
} %>% data.frame()

res_cmd <- rbind('conda activate mageckenv', res_cmd, 'conda deactivate')
write_tsv(res_cmd, file.path(out_dir, 'run_mageck_test.sh'), col_names = FALSE)
```

#### Load MAGeCK results

``` r
gene_summ_files <- list.files(path = file.path(params$output_dir, params$ne_mageck_dir), pattern = 'gene_summary.txt', full.names = TRUE)
analysis_prefix <- basename(gene_summ_files) %>% gsub(".gene_summary.txt", "", .)
ne_gene_summ <- foreach(i = gene_summ_files) %do% read.delim(i, check.names = FALSE)
names(ne_gene_summ) <- analysis_prefix

sgrna_summ_files <- list.files(path = file.path(params$output_dir, params$ne_mageck_dir), pattern = 'sgrna_summary.txt', full.names = TRUE)
analysis_prefix <- basename(sgrna_summ_files) %>% gsub(".sgrna_summary.txt", "", .)
ne_sgrna_summ <- foreach(i = sgrna_summ_files) %do% read.delim(i, check.names = FALSE)
names(ne_sgrna_summ) <- analysis_prefix
```


## sgRNA distribution

``` r
selected_panels <- list(
  MVT1 = c(
    rev(c('CRISPRKO_Clone_6', 'CRISPRKO_Clone_5', 'CRISPRKO_Clone_7', 'CRISPRKO_Clone_8', 'CRISPRKO_Clone_9')),
    rev(c('CRISPRKO_Clone_23', 'CRISPRKO_Clone_20', 'CRISPRKO_Clone_21', 'CRISPRKO_Clone_24', 'CRISPRKO_Clone_22'))
    ),
  `4T1` = c(
    rev(c('CRISPRKO_Clone_46', 'CRISPRKO_Clone_47', 'CRISPRKO_Clone_49', 'CRISPRKO_Clone_45', 'CRISPRKO_Clone_48')),
    rev(c('CRISPRKO_Clone_12', 'CRISPRKO_Clone_10', 'CRISPRKO_Clone_11', 'CRISPRKO_Clone_14', 'CRISPRKO_Clone_13'))
    ),
  `Py2T` = c(
    rev(c('CRISPRKO_Clone_50', 'CRISPRKO_Clone_51', 'CRISPRKO_Clone_52', 'CRISPRKO_Clone_53', 'CRISPRKO_Clone_54')),
    rev(c('CRISPRKO', 'CRISPRKO_Clone_1', 'CRISPRKO_Clone_2', 'CRISPRKO_Clone_3', 'CRISPRKO_Clone_4'))
    )
)

selected_panels_df <- foreach(i = names(selected_panels), .combine = rbind) %do% {
  data.frame(group = i, sample_alias = selected_panels[[i]])
}

use_samples <- unlist(selected_panels)
use_assay <- assay(se[,use_samples], 'cpm') %>% 
  as.data.frame(check.names = F) %>% 
  rownames_to_column('guide') %>% 
  pivot_longer(-guide, names_to = 'sample_alias', values_to = 'cpm') %>% 
  left_join(colData(se)%>% data.frame) %>% 
  left_join(selected_panels_df) %>% 
  mutate(
    sample_alias = factor(sample_alias, levels = use_samples),
    group = factor(group, levels = names(selected_panels))
  )
```

### Figure 4B: sgRNA distribution in mice transplanted with 4T1 cells previously subjected to StealTHY KO

Density plots showing unique sgRNAs from 4T1 orthotopic tumors (n=5-7); frequency (cpm) = read counts per million.


``` r
use_samples <- selected_panels_df %>% 
  filter(group == '4T1') %>% 
  pull(sample_alias)

use_assay <- assay(se[,use_samples], 'cpm') %>% 
  as.data.frame(check.names = F) %>% 
  rownames_to_column('guide') %>% 
  pivot_longer(-guide, names_to = 'sample_alias', values_to = 'cpm') %>% 
  left_join(colData(se)%>% data.frame) %>% 
  left_join(selected_panels_df) %>% 
  mutate(
    sample_alias = factor(sample_alias, levels = use_samples),
    group = factor(group, levels = names(selected_panels))
  )


use_assay %>% 
  ggplot(aes(x=cpm + 1, y = sample_alias, 
             fill = mouse_model_short, 
             color = mouse_model_short, 
             height = after_stat(density))) +
  geom_density_ridges(linewidth = default_linewidth, 
                      alpha = 0.8, 
                      scale = 4,
                      rel_min_height = 0.001 # set the `rel_min_height` argument to remove tails,
                      # stat = "density"
                      ) +
  scale_fill_manual(values = palette_mmodel) +
  scale_color_manual(values = palette_mmodel_line) +
  scale_x_log10(expand = expansion(mult = c(0, 0)),
                breaks = c(1, 1000, 100000)
                ) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.2))) +
  labs(y = '') +
  guides(fill = 'none', color = 'none') +
  theme_ridges(font_size = 8, grid = F) +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_line(linewidth = default_linewidth, color = 'black'),
    axis.ticks.x = element_line(linewidth = default_linewidth, color = 'black'),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=4),
    axis.title.x = element_text(size=4),
    strip.text.x = element_text(size=4, hjust = 0)
  )
```

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/test-dist-stealthy-ko-ridges_plot-4t1-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-test-dist-stealthy-ko-ridges_plot-4t1-1">
  Past versions of test-dist-stealthy-ko-ridges_plot-4t1-1.png
  </button>
  </p>

  <div id="fig-test-dist-stealthy-ko-ridges_plot-4t1-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/test-dist-stealthy-ko-ridges_plot-4t1-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

### Supplementary figure 4B: sgRNA distribution in mice transplanted with 4T1 cells previously subjected to StealTHY KO

Density plots showing the distribution of unique sgRNAs retrieved in primary tumors from the indicated mice transplanted with MVT1 cells (left) or Py2T cells (right) previously subjected to Stealthy KO; relative frequency is shown as read counts per million (cpm+1); n=5-7 each; representative tumors are shown.


``` r
use_samples <- selected_panels_df %>% 
  filter(group %in% c('MVT1', 'Py2T')) %>% 
  pull(sample_alias)

use_assay <- assay(se[,use_samples], 'cpm') %>% 
  as.data.frame(check.names = F) %>% 
  rownames_to_column('guide') %>% 
  pivot_longer(-guide, names_to = 'sample_alias', values_to = 'cpm') %>% 
  left_join(colData(se)%>% data.frame) %>% 
  left_join(selected_panels_df) %>% 
  mutate(
    sample_alias = factor(sample_alias, levels = use_samples),
    group = factor(group, levels = names(selected_panels))
  )


use_assay %>% 
  ggplot(aes(x=cpm + 1, y = sample_alias, 
             fill = mouse_model_short, 
             color = mouse_model_short, 
             height = after_stat(density))) +
  geom_density_ridges(linewidth = default_linewidth, 
                      alpha = 0.8, 
                      scale = 4,
                      rel_min_height = 0.001 # set the `rel_min_height` argument to remove tails,
                      ) +
  scale_fill_manual(values = palette_mmodel) +
  scale_color_manual(values = palette_mmodel_line) +
  scale_x_log10(expand = expansion(mult = c(0, 0)),
                breaks = c(1, 1000, 100000)
                ) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.2))) +
  labs(y = '') +
  guides(fill = 'none', color = 'none') +
  ggh4x::facet_wrap2(vars(group), nrow = 1,
                     scales = 'free_y', 
                     axes = "x") +
  theme_ridges(font_size = 8, grid = F) +
  theme(
    strip.background = element_blank(),
    axis.line.x = element_line(linewidth = default_linewidth, color = 'black'),
    axis.ticks.x = element_line(linewidth = default_linewidth, color = 'black'),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=4),
    axis.title.x = element_text(size=4),
    strip.text.x = element_text(size=4, hjust = 0)
  )
```

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/test-dist-stealthy-ko-ridges_plot-mvt1-py2t-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-test-dist-stealthy-ko-ridges_plot-mvt1-py2t-1">
  Past versions of test-dist-stealthy-ko-ridges_plot-mvt1-py2t-1.png
  </button>
  </p>

  <div id="fig-test-dist-stealthy-ko-ridges_plot-mvt1-py2t-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/test-dist-stealthy-ko-ridges_plot-mvt1-py2t-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
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
    scale_color_manual(values = palette_vector_mmodel) +
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

### Supplementary Figure 4D : K-S test table 

Table depicting the results of Kolmogorov-Smirnov tests, comparing sgRNA distributions from the CRISPR KO datasets indicated, encompassing three mammary carcinoma models (MVT1, 4T1, Py2T). The datasets are obtained from primary tumors harvested from the experiments shown in Figure 4A and Supplementary Figure 4A; n=5-7 each.


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
<div class="datatables html-widget html-fill-item" id="htmlwidget-832cd5f680003d8dd78c" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-832cd5f680003d8dd78c">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.21920240782541\" data-max=\"0.382583396037086\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>Results of Kolmogorov-Smirnov test<\/caption>","data":[["StealTHY syngeneic vs. StealTHY NSG","StealTHY syngeneic vs. Cas9+Puro"],[0.2192024078254109,0.3825833960370851],["&lt; 2.22e-16","&lt; 2.22e-16"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>comparison<\/th>\n      <th>D<\/th>\n      <th>P<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) {\n    return type !== 'display' ? data : DTWidget.formatRound(data, 3, 3, \",\", \".\", null);\n  }"},{"width":120,"targets":1},{"className":"dt-right","targets":1},{"name":"comparison","targets":0},{"name":"D","targets":1},{"name":"P","targets":2}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":["options.columnDefs.0.render"],"jsHooks":[]}</script>
```


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

### Dip test table for Syngeneic and NSG models

``` r
modality_test_df %>% 
  filter(!is.na(modality_group)) %>% 
  filter(modality_group != "StealTHY syngeneic") %>%
  mutate(
    modality_group = factor(
      modality_group, 
      levels = c("StealTHY syngeneic 4T1",
                 "Cas9+Puro syngeneic",
                 "StealTHY NSG",
                 "Cas9+Puro NSG")
      )
  ) %>%
  mutate(modality_group = fct_recode(
    modality_group, 
    `StealTHY syngeneic` = "StealTHY syngeneic 4T1")
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
  kbl(caption = 'Dip test results for Syngeneic and NSG models') %>%
  kable_paper(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)
```

<table class=" lightable-paper" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
<caption>Dip test results for Syngeneic and NSG models</caption>
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
   <td style="text-align:left;"> StealTHY syngeneic </td>
   <td style="text-align:right;"> 0.019 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 60.0 </td>
   <td style="text-align:right;"> 40.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cas9+Puro syngeneic </td>
   <td style="text-align:right;"> 0.068 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 100.0 </td>
   <td style="text-align:right;"> 0.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> StealTHY NSG </td>
   <td style="text-align:right;"> 0.031 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 86.7 </td>
   <td style="text-align:right;"> 13.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cas9+Puro NSG </td>
   <td style="text-align:right;"> 0.021 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 80.0 </td>
   <td style="text-align:right;"> 20.0 </td>
  </tr>
</tbody>
</table>


## Differential abundance analysis

``` r
gene_info <- readxl::read_xlsx('./data/resources/selected_gene_list/gene_list_for_figure_5c.xlsx')

use_colors <- c(
  `Immune system cytokines and chemokines` = 'grey',
  `Immune system cytokines\nand chemokines` = 'grey',
  
  `Interferone, TNF and danger/death-inducing signals` = 'black',
  `Interferon, TNF\nand danger/death-inducing\nsignals` = 'black',
  `Interferon, TNF and\ndanger/death-inducing\nsignals` = 'black',
  
  `TGF-β family related` = 'red',
  
  `RTK signaling` = 'blue',
  `Positive controls` = 'black'
)

use_fill <- c(
  `Immune system cytokines and chemokines` = 'grey',
  `Immune system cytokines\nand chemokines` = 'grey',
  
  `Interferone, TNF and danger/death-inducing signals` = 'black',
  `Interferon, TNF\nand danger/death-inducing\nsignals` = 'black',
  `Interferon, TNF and\ndanger/death-inducing\nsignals` = 'black',
  `TGF-β family related` = 'red',
  
  `RTK signaling` = 'blue',
  `Positive controls` = 'white'
)

use_shape <- c(
  `Immune system cytokines and chemokines` = 16,
  `Immune system cytokines\nand chemokines` = 16,
  
  `Interferone, TNF and danger/death-inducing signals` = 16,
  `Interferon, TNF\nand danger/death-inducing\nsignals` = 16,
  `Interferon, TNF and\ndanger/death-inducing\nsignals` = 16,
  `TGF-β family related` = 16,
  
  `RTK signaling` = 16,
  `Positive controls` = 1
)

# get gene order based on Comparison : `StealTHY KO Syngeneic` 
gene_levels <- gene_summ[['stealthy_ko_syngeneic']] %>%
  ReadRRA(score = 'lfc') %>%
  dplyr::rename(LFC = Score) %>%
  left_join(gene_info, by = c('id' = 'id')) %>%
  filter(!is.na(gene_set)) %>%
  mutate(gene_name = fct_reorder(gene_name, LFC, .desc = FALSE)) %>%
  pull(gene_name) %>%
  levels
```

### Figure 4C and Supplementary Figure 4F

Dot plots highlighting hits identified following the strategy depicted in Figure 4A upon transplantation in the indicated hosts of metastatic mammary carcinoma cells (MVT1, 4T1, Py2T) previously subjected to Stealthy KO or gold standard KO by Cas9+PuroR; representative genes from four recurrently identified gene families are shown in comparison across screening methods. The fold-change values indicate differential sgRNA abundance between primary tumors and synchronous lung metastases; n=5-7 each; RTK, receptor tyrosine kinase, Tnf, tumor necrosis factor; Ifn, interferon; FDR, false discovery rate.


``` r
use_comp <- c('stealthy_ko_syngeneic', 'stealthy_ko_nsg', 'cas9+puro_syngeneic')

# FDR limits
xfdr <- foreach(i=gene_summ[use_comp], .combine = rbind) %do% {
  i %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info, by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    pull(FDR) %>% 
    min
}
maxFDR <- -log10(xfdr[,1]) %>% max
 
# Set LFC (x-axis) for comparable plots to "stealthy_ko_syngeneic": "stealthy_ko_nsg"
stealthy_comp <- c("stealthy_ko_syngeneic", "stealthy_ko_nsg")
xlfc <- foreach(i=gene_summ[stealthy_comp], .combine = rbind) %do% {
  i %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info, by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    pull(LFC) %>% 
    abs %>% max
}
stealthy_maxLFC <- max(xlfc[,1])
stealthy_xlim <- c(-stealthy_maxLFC, stealthy_maxLFC)
stealthy_breaks <- c(-2.5, -1, 0, 1, 2.5)
stealthy_labels <- c('-2.5', '-1', '0', '1', '2.5')

i <- use_comp[1]
for(i in use_comp)  {
  cat("####", mageck_comp_list_names[i], "\n\n")
  
  x <- gene_summ[[i]] %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info, by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    mutate(
      gene_name = factor(gene_name, levels = gene_levels),
      gene_set = ifelse(
        gene_set == 'Immune system cytokines and chemokines',
        'Immune system cytokines\nand chemokines',
        gene_set
      ),
      gene_set = ifelse(
        gene_set == 'Interferone, TNF and danger/death-inducing signals',
        'Interferon, TNF and\ndanger/death-inducing\nsignals',
        gene_set
      ),
      gene_set = factor(gene_set, names(use_colors))
      )
  xlim <- c(-ceiling(max(abs(x$LFC), na.rm = TRUE)), 
            ceiling(max(abs(x$LFC), na.rm = TRUE)))
  size_breaks <- c(0, 0.25, 0.5, 1, 2, ceiling(maxFDR))
  res <- x %>% 
    ggplot(aes(LFC, gene_name, color = gene_set, size = -log10(FDR))) +
    geom_point(aes(shape = gene_set)) +
    geom_vline(xintercept = 0, linetype = 3, linewidth = one_pt/2) +
    scale_y_discrete(expand = expansion(add = c(1, 1))) + # add fixed units around each facet to avoid lines over points
    scale_color_manual(values = use_colors) +
    scale_size(range = c(0.5, 2.5), 
               limits = c(0, ceiling(maxFDR)), 
               breaks = size_breaks,
               labels = size_breaks %>% as.character # remove trailing 0 from decimals
               ) +
    scale_shape_manual(values = use_shape) +
    scale_x_continuous(
        expand = expansion(add = c(1, 1))
        ) +
    facet_grid(rows = vars(gene_set), scales = 'free_y', space = 'free_y') +
    theme_facet +
    theme(
      axis.ticks = element_line(linewidth = one_pt/2, color = 'black'),
      panel.border = element_rect(linewidth = one_pt/2, color = 'black'),
      axis.text = element_text(size=8), 
      axis.title = element_text(size = 8),
      legend.title = element_text(size=8),
      legend.text = element_text(size=8),
      strip.text = element_text(size = 4)
    ) +
    labs(
      x = expression(paste("lo", g[2],"(Fold change)")),
      y = NULL,
      size = expression(paste("-lo", g[10],"(FDR)")),
      title = 
    ) +
    guides(color = 'none', shape = 'none')
  
  if(i %in% stealthy_comp)
    res <- res +
      scale_x_continuous(
        expand = expansion(add = c(1, 1)), 
        limits = stealthy_xlim,
        breaks = stealthy_breaks,
        labels = stealthy_labels
      )
  
  print(res)
  cat('\n\n')
}
```

#### StealTHY KO Syngeneic 

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-da-selected-genes-dotplot-1">
  Past versions of da-selected-genes-dotplot-1.png
  </button>
  </p>

  <div id="fig-da-selected-genes-dotplot-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

#### StealTHY KO NSG 

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-2.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-da-selected-genes-dotplot-2">
  Past versions of da-selected-genes-dotplot-2.png
  </button>
  </p>

  <div id="fig-da-selected-genes-dotplot-2" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-2.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

#### Cas9+Puro Syngeneic 

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-3.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-da-selected-genes-dotplot-3">
  Past versions of da-selected-genes-dotplot-3.png
  </button>
  </p>

  <div id="fig-da-selected-genes-dotplot-3" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-selected-genes-dotplot-3.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


### Supplementary Figure 4E

Dot plots presenting differential sgRNA abundance for selected positive control genes in mammary carcinoma models (MVT1, 4T1, Py2T) subjected to StealTHY KO with respect to a mock KO reference dataset. In the mock KO condition, dCas9 transduction was performed instead of Apo-Cas9 transfection in the same sgRNA-transduced cells. Top: dot plot showing differential sgRNA abundance measured in vitro in cells prior to transplantation, i.e. harvested 60 hours after Apo-Cas9 transfection or dCas9 transduction; n=3 samples per group. Bottom: dot plot showing differential sgRNA abundance measured in primary tumors in vivo, harvested at the experimental endpoint from syngeneic hosts; n=5-7 mice per group.


``` r
use_comp <- c('stealthy_ko_over_dcas9_in_vitro', 'stealthy_ko_over_dcas9_in_vivo')

# FDR limits
xfdr <- foreach(i=gene_summ[use_comp], .combine = rbind) %do% {
  i %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info, by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    pull(FDR) %>% 
    min
}
maxFDR <- -log10(xfdr[,1]) %>% max

i <- use_comp[1]
for(i in use_comp)  {
  cat("####", mageck_comp_list_names[i], "\n\n")
  
  x <- gene_summ[[i]] %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info %>% filter(gene_set == 'Positive controls'), 
              by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    mutate(
      gene_name = factor(gene_name, levels = gene_levels),
      gene_set = ifelse(
        gene_set == 'Immune system cytokines and chemokines',
        'Immune system cytokines\nand chemokines',
        gene_set
      ),
      gene_set = ifelse(
        gene_set == 'Interferone, TNF and danger/death-inducing signals',
        'Interferon, TNF and\ndanger/death-inducing\nsignals',
        gene_set
      ),
      gene_set = factor(gene_set, names(use_colors))
      )
  xlim <- c(-ceiling(max(abs(x$LFC), na.rm = TRUE)), 
            ceiling(max(abs(x$LFC), na.rm = TRUE)))
  size_breaks <- c(0, 0.25, 0.5, 1, 2, ceiling(maxFDR))
  res <- x %>% 
    ggplot(aes(LFC, gene_name, size = -log10(FDR))) +
    geom_point(shape = 1, linewidth = default_linewidth) +
    
    geom_vline(xintercept = 0, linetype = 3, linewidth = one_pt/2) +
    scale_y_discrete(expand = expansion(add = c(1, 1))) + # add fixed units around each facet to avoid lines over points
    scale_color_manual(values = use_colors) +
    scale_size(range = c(2, 4), 
               limits = c(0, ceiling(maxFDR)), 
               breaks = size_breaks,
               labels = size_breaks %>% as.character # remove trailing 0 from decimals
               ) +
    scale_shape_manual(values = use_shape) +
    scale_x_continuous(
        expand = expansion(add = c(1, 1)),
        limits = c(-4.5, 4.5),
        breaks = c(-4, 0, 4),
        ) +
    facet_grid(rows = vars(gene_set), scales = 'free_y', space = 'free_y') +
    theme_facet +
    theme(
      axis.ticks = element_line(linewidth = one_pt/2, color = 'black'),
      panel.border = element_rect(linewidth = one_pt/2, color = 'black'),
      axis.text = element_text(size=8), 
      axis.title = element_text(size = 8),
      legend.title = element_text(size=8),
      legend.text = element_text(size=8),
      strip.text = element_text(size = 4)
    ) +
    labs(
      x = expression(paste("lo", g[2],"(Fold change)")),
      y = NULL,
      size = expression(paste("-lo", g[10],"(FDR)")),
      title = 
    ) +
    guides(color = 'none', shape = 'none')
  
  print(res)
  cat('\n\n')
}
```

#### StealTHY KO over dCas9 In vitro 

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-positive-controls-genes-dotplot-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-da-positive-controls-genes-dotplot-1">
  Past versions of da-positive-controls-genes-dotplot-1.png
  </button>
  </p>

  <div id="fig-da-positive-controls-genes-dotplot-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-positive-controls-genes-dotplot-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

#### StealTHY KO over dCas9 In vivo 

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-positive-controls-genes-dotplot-2.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-da-positive-controls-genes-dotplot-2">
  Past versions of da-positive-controls-genes-dotplot-2.png
  </button>
  </p>

  <div id="fig-da-positive-controls-genes-dotplot-2" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/da-positive-controls-genes-dotplot-2.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  



## Non-expressed hits

- Definition of **non-expressed genes** : We keep genes with ≥ 1 read count in 1/3 of the samples for each site (whole tumor or whole lung). In this case, if only 3 samples for each site are available, the gene must have at least 1 read count in 1 sample. Then, the final list of non-expressed genes is defined as the intersection of non-expressed genes in whole tumor and whole lung.

- Definition of **hits** : Hits are defined as genes with a log2 fold-change ≤ -2 or ≥ 2.


``` r
# Color palette
use_color_palette <- c(
  `StealTHY`  = '#F7EC13',
  `Puro` = 'black',
  `Cas9+Puro` = 'black',
  `StealTHY NSG`  = '#58595B'
)
```


``` r
min_counts <- 1
min_prop <- 1/3

# Read RNA-seq data
se_rna <- readRDS(params$se_rnaseq_facs)

# Calculate CPM
assay(se_rna, 'cpm') <- edgeR::cpm(se_rna)

# CRISPR library genes
crispr_genes <- unique(rowData(se)$Gene) %>% grep('Non_Target', ., invert = T, value = T)

# Non-expressed genes
non_expr_genes <- list()
expr_genes <- list()
```


``` r
use_samples <- c('4T1_In_vivo_Thy1', 'MVT1_In_vivo_Thy1', 'Py2T_In_vivo_Thy1')

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <- rowData(use_se)$gene_name

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_tumor']] <- rows_to_keep
non_expr_genes[['whole_tumor']] <- rows_to_remove
```


``` r
use_samples <- c('4T1_In_vivo_Thy1_Lungs', 'MVT1_In_vivo_Thy1_Lungs', 'Py2T_In_vivo_Thy1_Lungs')

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <- rowData(use_se)$gene_name

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_lung']] <- rows_to_keep
non_expr_genes[['whole_lung']] <- rows_to_remove
```



``` r
use_annot <- colData(se)

expr_genes_combined <- expr_genes
non_expr_genes_combined <- non_expr_genes

use_comp <- c("StealTHY", "Puro", "Cas9+Puro", "StealTHY NSG")
i <- 'StealTHY'
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes[[st1]], expr_genes[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_global <- hits_stats
```



``` r
min_counts <- 1
min_prop <- 1/3
use_donor <- '4T1'
whole_tumor_sample <- '4T1_In_vivo_Thy1'
whole_lung_sample <- '4T1_In_vivo_Thy1_Lungs'

# Read RNA-seq data
se_rna <- readRDS(params$se_rnaseq_facs)

# Calculate CPM
assay(se_rna, 'cpm') <- edgeR::cpm(se_rna)

# CRISPR library genes
crispr_genes <- unique(rowData(se)$Gene) %>% grep('Non_Target', ., invert = T, value = T)

# Non-expressed genes
non_expr_genes <- list()
expr_genes <- list()
```


``` r
use_samples <- whole_tumor_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_tumor']] <- rows_to_keep
non_expr_genes[['whole_tumor']] <- rows_to_remove
```


``` r
use_samples <- whole_lung_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_lung']] <- rows_to_keep
non_expr_genes[['whole_lung']] <- rows_to_remove
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY 4T1", "Puro 4T1", "Cas9+Puro 4T1", "StealTHY NSG 4T1")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes[[st1]], expr_genes[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / prop_non_expr_genes,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_4t1 <- hits_stats %>% 
  mutate(
    donor = '4T1',
    comparison = gsub(" 4T1", "", comparison)
    )
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY 4T1", "Puro 4T1", "Cas9+Puro 4T1", "StealTHY NSG 4T1")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes_combined[[st1]], expr_genes_combined[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / prop_non_expr_genes,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_4t1_comb <- hits_stats %>% 
  mutate(
    donor = '4T1',
    comparison = gsub(" 4T1", "", comparison)
    )
```




``` r
min_counts <- 1
min_prop <- 1/3
use_donor <- 'mvt1'
whole_tumor_sample <- 'MVT1_In_vivo_Thy1'
whole_lung_sample <- 'MVT1_In_vivo_Thy1_Lungs'

# Read RNA-seq data
se_rna <- readRDS(params$se_rnaseq_facs)

# Calculate CPM
assay(se_rna, 'cpm') <- edgeR::cpm(se_rna)

# CRISPR library genes
crispr_genes <- unique(rowData(se)$Gene) %>% grep('Non_Target', ., invert = T, value = T)

# Non-expressed genes
non_expr_genes <- list()
expr_genes <- list()
```


``` r
use_samples <- whole_tumor_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_tumor']] <- rows_to_keep
non_expr_genes[['whole_tumor']] <- rows_to_remove
```


``` r
use_samples <- whole_lung_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_lung']] <- rows_to_keep
non_expr_genes[['whole_lung']] <- rows_to_remove
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY MVT1", "Puro MVT1", "StealTHY NSG MVT1")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes[[st1]], expr_genes[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / n_nohits,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_mvt1 <- hits_stats %>% 
  mutate(
    donor = 'MVT1',
    comparison = gsub(" MVT1", "", comparison)
    )
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY MVT1", "Puro MVT1", "StealTHY NSG MVT1")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes_combined[[st1]], expr_genes_combined[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / n_nohits,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_mvt1_comb <- hits_stats %>% 
  mutate(
    donor = 'MVT1',
    comparison = gsub(" MVT1", "", comparison)
    )
```



``` r
min_counts <- 1
min_prop <- 1/3
use_donor <- 'py2t'
whole_tumor_sample <- 'Py2T_In_vivo_Thy1'
whole_lung_sample <- 'Py2T_In_vivo_Thy1_Lungs'

# Read RNA-seq data
se_rna <- readRDS(params$se_rnaseq_facs)

# Calculate CPM
assay(se_rna, 'cpm') <- edgeR::cpm(se_rna)

# CRISPR library genes
crispr_genes <- unique(rowData(se)$Gene) %>% grep('Non_Target', ., invert = T, value = T)

# Non-expressed genes
non_expr_genes <- list()
expr_genes <- list()
```


``` r
use_samples <- whole_tumor_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_tumor']] <- rows_to_keep
non_expr_genes[['whole_tumor']] <- rows_to_remove
```


``` r
use_samples <- whole_lung_sample

# RNA-seq : Select genes and samples
use_se <- se_rna[rowData(se_rna)$gene_name %in% crispr_genes, use_samples]

# Change rownames
rownames(use_se) <-make.names(rowData(use_se)$gene_name, unique = TRUE)

# Proportion of samples with min_counts at each group
use_group <- as.factor(use_se$donor)
counts_mat <- assay(use_se, 'counts')
counts_stats <- foreach(i=levels(use_group), .combine = cbind) %do% {
  use_cols <- use_se$donor == i
  rowSums(counts_mat[,use_cols, drop = FALSE] >= min_counts)  / sum(use_cols)
} %>% data.frame
colnames(counts_stats) <- levels(use_group)
counts_stats$all <- rowSums(counts_mat >= min_counts)  / ncol(counts_mat)
  
# Keep genes with min_counts threshold in >=min_prop of samples in at least 1 group
rows_to_keep <- rownames(counts_stats)[counts_stats$all >= min_prop]
rows_to_remove <- setdiff(rownames(counts_stats), rows_to_keep)

expr_genes[['whole_lung']] <- rows_to_keep
non_expr_genes[['whole_lung']] <- rows_to_remove
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY Py2T", "Puro Py2T", "StealTHY NSG Py2T")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes[[st1]], expr_genes[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / n_nohits,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_py2t <- hits_stats %>% 
  mutate(
    donor = 'Py2T',
    comparison = gsub(" Py2T", "", comparison)
  )
```


``` r
use_annot <- colData(se)

use_comp <- c("StealTHY Py2T", "Puro Py2T", "StealTHY NSG Py2T")

i <- use_comp[1]
hits_stats <- foreach(i = use_comp, .combine = rbind) %do% {
  i_names <- i %>% gsub(" ", "_", .) %>% tolower
  st1 <- unique(use_annot[ne_mageck_comp_list[[i]][[1]],]$sample_type)
  st2 <- unique(use_annot[ne_mageck_comp_list[[i]][[2]],]$sample_type)
  
  genes_expr <- intersect(expr_genes_combined[[st1]], expr_genes_combined[[st2]])
 
  da_res <- ne_gene_summ[[i_names]]  %>%  
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    mutate(
      expressed = id %in% genes_expr,
      hit = (abs(LFC) >= 2)
      )
  
   prop_non_expr_genes <- (nrow(da_res)-length(genes_expr)) / nrow(da_res)
   prop_non_expr_hits <- sum(da_res$hit & !da_res$expressed) / sum(da_res$hit)
   data.frame(
     comparison = i,
     total_genes = nrow(da_res),
     n_nonexpr_genes = nrow(da_res)-length(genes_expr),
     prop_non_expr_genes = (nrow(da_res)-length(genes_expr)) / nrow(da_res),
     n_hits = sum(da_res$hit),
     non_expr_hits = sum(da_res$hit & !da_res$expressed),
     n_nohits = sum(!da_res$hit)
   ) %>% 
     mutate(
       prop_non_expr_hits = non_expr_hits / n_hits,
       prop_non_expr_hits_over_non_expr = non_expr_hits / n_nohits,
       lfc = log2(prop_non_expr_hits / prop_non_expr_genes)
     )
  
}

hits_stats <- hits_stats %>% 
  mutate(
    comparison = factor(comparison,levels = names(ne_mageck_comp_list) %>% rev)
  )

hits_stats_py2t_comb <- hits_stats %>% 
  mutate(
    donor = 'Py2T',
    comparison = gsub(" Py2T", "", comparison)
  )
```



### Figure 4D: Proportion of non-expressed hits

Bar graph presenting the proportion of hits (log2 fold-change ≥ 2 or ≤ −2) not expressed (read counts <1 in ≥ 2/3 of samples) in MVT1, 4T1 and Py2T models (n=5-7), assessed using RNA-seq of Thy(A)-sorted carcinoma cells from primary tumors (n=3) and metastatic lungs (n=3).


``` r
chunk_data <- hits_stats_global

chunk_data %>%
  ggplot(aes(comparison, prop_non_expr_hits, fill = comparison)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 'solid',
             linewidth = one_pt/4) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
    )  +
  scale_fill_manual(values = use_color_palette) +
  theme_facet +
  theme(
    axis.line = element_line(linewidth = one_pt/4, color = 'black'),
    axis.ticks = element_line(linewidth = one_pt/4, color = 'black'),
    panel.border = element_blank(),
    axis.text = element_text(size=5), 
    axis.title = element_text(size = 5),
    legend.title = element_text(size=5),
    legend.text = element_text(size=5),
    strip.text.y = element_text(size = 4, angle = 0)
  ) +
  guides(fill = 'none') +
  labs(
    x = '',
    y = 'Proportion of\nnon-expressed hits'
  )
```

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/non-expressed-hits-barplot-prop-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-non-expressed-hits-barplot-prop-1">
  Past versions of non-expressed-hits-barplot-prop-1.png
  </button>
  </p>

  <div id="fig-non-expressed-hits-barplot-prop-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/non-expressed-hits-barplot-prop-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


### Supplementary figure 4G: Proportion of non-expressed hits

Detailed analysis of non-expressed screen hits sub-divided per individual mammary carcinoma model. Top panel: bar graph presenting the proportion of hits that are not expressed but are still identified by CRISPR KO screening. Genes were classified as “hit” if their absolute log2 fold-change was ≥2 in each screening condition (minimum n=5 mice per model), and as “non-expressed” if read counts were <1 in at least 2/3 of the samples, according to RNA-seq of Thy(A)-sorted carcinoma cells from metastatic lungs (n=3) in intersection with RNA-seq data from primary tumors (n=3); all RNA-seq data have been acquired in syngeneic hosts.


``` r
chunk_data <- rbind(hits_stats_4t1, hits_stats_mvt1, hits_stats_py2t) %>% 
  mutate(
    comparison = factor(comparison, levels = rev(c('StealTHY', 'Puro', 'Cas9+Puro', 'StealTHY NSG')))
  )
chunk_data %>%
  ggplot(aes(comparison, prop_non_expr_hits, fill = comparison)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
    ) +
  scale_fill_manual(values = use_color_palette) +
  facet_grid(rows = vars(donor),
             scales = 'free_y',
             space = 'free_y') +
  theme_facet +
  theme(
    axis.ticks = element_line(linewidth = one_pt/2, color = 'black'),
    panel.border = element_rect(linewidth = one_pt/2, color = 'black'),
    axis.text = element_text(size=4), 
    axis.title = element_text(size = 4),
    legend.title = element_text(size=4),
    legend.text = element_text(size=4),
    strip.text.y = element_text(size = 4, angle = 0)
  ) +
  guides(fill = 'none') +
  labs(
    x = '',
    y =  'Proportion of\nnon-expressed hits'
  )
```

<img src="figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/non-expressed-hits-barplot-prop-bymodel-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-non-expressed-hits-barplot-prop-bymodel-1">
  Past versions of non-expressed-hits-barplot-prop-bymodel-1.png
  </button>
  </p>

  <div id="fig-non-expressed-hits-barplot-prop-bymodel-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/32aca023d6d8209f1aa3df3f177fc383c7569844/docs/figure/crispr-mm_2215_sgRNA-StealTHY.Rmd/non-expressed-hits-barplot-prop-bymodel-1.png" target="_blank">32aca02</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

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
 [1] cowplot_1.2.0               MAGeCKFlute_2.9.0          
 [3] ggpubr_0.6.1                ggbeeswarm_0.7.2           
 [5] colorblindr_0.1.0           colorspace_2.1-1           
 [7] patchwork_1.3.1             ggh4x_0.3.1                
 [9] ggridges_0.5.6              SummarizedExperiment_1.36.0
[11] Biobase_2.66.0              GenomicRanges_1.58.0       
[13] GenomeInfoDb_1.42.3         IRanges_2.40.1             
[15] S4Vectors_0.44.0            BiocGenerics_0.52.0        
[17] MatrixGenerics_1.18.1       matrixStats_1.5.0          
[19] diptest_0.77-1              kableExtra_1.4.0           
[21] DT_0.34.0                   magrittr_2.0.3             
[23] foreach_1.5.2               knitr_1.50                 
[25] lubridate_1.9.4             forcats_1.0.0              
[27] stringr_1.5.1               dplyr_1.1.4                
[29] purrr_1.1.0                 readr_2.1.5                
[31] tidyr_1.3.1                 tibble_3.3.0               
[33] ggplot2_3.5.2               tidyverse_2.0.0            
[35] workflowr_1.7.1            

loaded via a namespace (and not attached):
  [1] splines_4.4.3           later_1.4.4             ggplotify_0.1.2        
  [4] bitops_1.0-9            filelock_1.0.3          cellranger_1.1.0       
  [7] R.oo_1.27.1             graph_1.84.1            XML_3.99-0.19          
 [10] lifecycle_1.0.4         httr2_1.2.1             rstatix_0.7.2          
 [13] edgeR_4.4.2             rprojroot_2.1.0         processx_3.8.6         
 [16] lattice_0.22-7          crosstalk_1.2.2         backports_1.5.0        
 [19] limma_3.62.2            sass_0.4.10             rmarkdown_2.29         
 [22] jquerylib_0.1.4         yaml_2.3.10             ggtangle_0.0.7         
 [25] httpuv_1.6.16           depmap_1.20.0           DBI_1.2.3              
 [28] RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.52.0        
 [31] R.utils_2.13.0          RCurl_1.98-1.17         yulab.utils_0.2.0      
 [34] rappdirs_0.3.3          git2r_0.36.2            GenomeInfoDbData_1.2.13
 [37] enrichplot_1.26.6       ggrepel_0.9.6           tidytree_0.4.6         
 [40] svglite_2.2.1           codetools_0.2-20        DelayedArray_0.32.0    
 [43] DOSE_4.0.1              xml2_1.3.8              tidyselect_1.2.1       
 [46] aplot_0.2.8             UCSC.utils_1.2.0        farver_2.1.2           
 [49] BiocFileCache_2.14.0    pathview_1.46.0         jsonlite_2.0.0         
 [52] Formula_1.2-5           iterators_1.0.14        systemfonts_1.2.3      
 [55] tools_4.4.3             treeio_1.30.0           Rcpp_1.1.0             
 [58] glue_1.8.0              gridExtra_2.3           SparseArray_1.6.2      
 [61] xfun_0.53               qvalue_2.38.0           withr_3.0.2            
 [64] BiocManager_1.30.26     fastmap_1.2.0           callr_3.7.6            
 [67] digest_0.6.37           gridGraphics_0.5-1      timechange_0.3.0       
 [70] R6_2.6.1                textshaping_1.0.1       GO.db_3.20.0           
 [73] dichromat_2.0-0.1       RSQLite_2.4.3           R.methodsS3_1.8.2      
 [76] generics_0.1.4          data.table_1.17.8       httr_1.4.7             
 [79] htmlwidgets_1.6.4       S4Arrays_1.6.0          whisker_0.4.1          
 [82] pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4             
 [85] XVector_0.46.0          clusterProfiler_4.14.6  htmltools_0.5.8.1      
 [88] carData_3.0-5           fgsea_1.32.4            scales_1.4.0           
 [91] png_0.1-8               ggfun_0.2.0             rstudioapi_0.17.1      
 [94] tzdb_0.5.0              reshape2_1.4.4          nlme_3.1-168           
 [97] curl_7.0.0              org.Hs.eg.db_3.20.0     cachem_1.1.0           
[100] BiocVersion_3.20.0      parallel_4.4.3          vipor_0.4.7            
[103] AnnotationDbi_1.68.0    pillar_1.11.0           grid_4.4.3             
[106] vctrs_0.6.5             promises_1.3.3          car_3.1-3              
[109] dbplyr_2.5.0            beeswarm_0.4.0          Rgraphviz_2.50.0       
[112] evaluate_1.0.5          KEGGgraph_1.66.0        locfit_1.5-9.12        
[115] cli_3.6.5               compiler_4.4.3          rlang_1.1.6            
[118] crayon_1.5.3            ggsignif_0.6.4          labeling_0.4.3         
[121] ps_1.9.1                getPass_0.2-4           plyr_1.8.9             
[124] fs_1.6.6                stringi_1.8.7           viridisLite_0.4.2      
[127] BiocParallel_1.40.2     Biostrings_2.74.1       lazyeval_0.2.2         
[130] GOSemSim_2.32.0         Matrix_1.7-3            ExperimentHub_2.14.0   
[133] hms_1.1.3               bit64_4.6.0-1           statmod_1.5.0          
[136] KEGGREST_1.46.0         AnnotationHub_3.14.0    igraph_2.1.4           
[139] broom_1.0.9             memoise_2.0.1           bslib_0.9.0            
[142] ggtree_3.14.0           fastmatch_1.1-6         bit_4.6.0              
[145] readxl_1.4.5            gson_0.1.0              ape_5.8-1              
</div>
