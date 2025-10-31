---
title: "StealTHY CRISPR screen"
subtitle: "Mouse libray (n = 2215) - Reviewer 2 Observation C"
author: "Francesc Castro-Giner"
date: 'October 31, 2025'
output: workflowr::wflow_html
editor_options:
  chunk_output_type: console
params:
  data_dir: ./data/crispr/mm_2215_sgRNA
  output_dir: ./output/crispr/mm_2215_sgRNA
  mageck_dir: mageck_rra_r2oC
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
	Untracked:  analysis/crispr-mm_2215_sgRNA-StealTHY.md
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
to the R Markdown (<code>analysis/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd</code>) and HTML (<code>docs/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.html</code>)
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
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/455fa9121fcca51363cf639f7ee170fb108e34c1/analysis/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd" target="_blank">455fa91</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Prepare files for publication</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.html" target="_blank">a6972e9</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/b6bc77088b14bf2712a35dc41be83f6032057e27/analysis/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd" target="_blank">b6bc770</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-07-01</td>
<td>update crispr analysis</td>
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
library(openxlsx)

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

### MAGeCK analysis

#### List of comparisons
List of comparisons for differential abundance using MAGeCK, using Mets over Primary tumor

``` r
x <- colData(se) %>% data.frame
mageck_comp_list <- list(
  `StealTHY KO Syngeneic CMT167` = list(
    `Liver mets` = x %>% filter(vector_mmodel == 'CMT167-KO-Thy1-C57BL6' & sample_type == 'Liver mets' ) %>% pull(sample_alias),
    `Primary tumor` = x %>% filter(vector_mmodel == 'CMT167-KO-Thy1-C57BL6' & sample_type == 'Primary tumor' ) %>% pull(sample_alias)
    ),
  `StealTHY KO Syngeneic LLC1` = list(
    `Liver mets` = x %>% filter(vector_mmodel == 'LLC1-KO-Thy1-C57BL6' & sample_type == 'Liver mets' ) %>% pull(sample_alias),
    `Primary tumor` = x %>% filter(vector_mmodel == 'LLC1-KO-Thy1-C57BL6' & sample_type == 'Primary tumor' ) %>% pull(sample_alias)
    ),
  `StealTHY KO Syngeneic CT26` = list(
    `Liver mets` = x %>% filter(vector_mmodel == 'CT26-KO-Thy1-BALB' & sample_type == 'Lung mets' ) %>% pull(sample_alias),
    `Primary tumor` = x %>% filter(vector_mmodel == 'CT26-KO-Thy1-BALB' & sample_type == 'Primary tumor' ) %>% pull(sample_alias)
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
<script type="application/json" data-for="htmlwidget-0e41329c1963262c826a">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"caption":"<caption>List of comparisons for MAGeCK analysis<\/caption>","data":[["StealTHY KO Syngeneic CMT167","StealTHY KO Syngeneic LLC1","StealTHY KO Syngeneic CT26"],["C167_StealTHY_Met_1 C167_StealTHY_Met_3 C167_StealTHY_Met_4 C167_StealTHY_Met_5","LLC1_StealTHY_Met_1 LLC1_StealTHY_Met_2 LLC1_StealTHY_Met_3 LLC1_StealTHY_Met_4 LLC1_StealTHY_Met_5","ColoCT26_StealTHY_Met_1 ColoCT26_StealTHY_Met_2 Replacement_1 Replacement_3 Replacement_4bis"],["C167_StealTHY_Tum_1 C167_StealTHY_Tum_2 C167_StealTHY_Tum_3 C167_StealTHY_Tum_4 C167_StealTHY_Tum_5","LLC1_StealTHY_Tum_1 LLC1_StealTHY_Tum_2 LLC1_StealTHY_Tum_3 LLC1_StealTHY_Tum_4 LLC1_StealTHY_Tum_5","ColoCT26_StealTHY_Tum_1 ColoCT26_StealTHY_Tum_2 ColoCT26_StealTHY_Tum_3 Replacement_2 Replacement_4"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Comparison<\/th>\n      <th>Case<\/th>\n      <th>Control<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["csv","excel"],"columnDefs":[{"width":120,"targets":[1,2]},{"name":"Comparison","targets":0},{"name":"Case","targets":1},{"name":"Control","targets":2}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>
```

#### Run MAGeCK

Generate MAGeCK test scripts. To generate the script files change `eval = FALSE` to `eval = TRUE` in the chunk options.

The scripts must be run in the terminal, inside the path . The script will activate the mageck environment and run the test for each comparison. The results will be saved in the same directory.

``` r
chunk_se <- se
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

gene_summ_files <- list.files(path = file.path(params$output_dir, params$mageck_dir), pattern = 'gene_summary.txt', full.names = TRUE)
analysis_prefix <- basename(gene_summ_files) %>% gsub(".gene_summary.txt", "", .)
gene_summ <- foreach(i = gene_summ_files) %do% read.delim(i, check.names = FALSE)
names(gene_summ) <- analysis_prefix


gene_summ_files <- list.files(path = file.path('./output/crispr/mm_2215_sgRNA', 'mageck_rra_expr'), pattern = 'gene_summary.txt', full.names = TRUE)
analysis_prefix <- basename(gene_summ_files) %>% gsub(".gene_summary.txt", "", .)
gene_summ_5C <- foreach(i = gene_summ_files) %do% read.delim(i, check.names = FALSE)
names(gene_summ_5C) <- analysis_prefix
```



## Download MAGeCK results

``` r
# File name gene summary
rmd_file <- current_input()
if(is.null(rmd_file))
  rmd_file <- 'tmp'
gene_file_xlsx <- file.path('./docs/file',rmd_file, 'mageck_gene_summary.xlsx')
sgrna_file_xlsx <- file.path('./docs/file',rmd_file, 'mageck_sgrna_summary.xlsx')
```


``` r
for(i in names(mageck_comp_list)) {
  i <- tolower(i) %>% gsub(" ", "_",.) 
  file_xlsx <- file.path('./docs/file',rmd_file, paste0(i, '_mageck_results.xlsx'))
  dir.create(dirname(file_xlsx), recursive = TRUE, showWarnings = FALSE)
  # Generate workbook
  wb <- createWorkbook()
  wsn <- 'gene_summ'
  addWorksheet(wb, wsn)
  writeData(wb, wsn, gene_summ[[i]])
  wsn <- 'sgrna_summ'
  addWorksheet(wb, wsn)
  saveWorkbook(wb, file_xlsx, TRUE)
}
```


``` r
# File name gene summary
file_xlsx <- gene_file_xlsx
dir.create(dirname(file_xlsx), recursive = TRUE, showWarnings = FALSE)

# Generate workbook
wb <- createWorkbook()
for(i in names(gene_summ)) {
  addWorksheet(wb, i)
  writeData(wb, i, gene_summ[[i]])
}
saveWorkbook(wb, file_xlsx, TRUE)

# File name sgRNA summary
file_xlsx <- sgrna_file_xlsx
dir.create(dirname(file_xlsx), recursive = TRUE, showWarnings = FALSE)

# Generate workbook
wb <- createWorkbook()
for(i in names(sgrna_summ)) {
  addWorksheet(wb, i)
  writeData(wb, i, sgrna_summ[[i]])
}
saveWorkbook(wb, file_xlsx, TRUE)
```

### Overall summaries
The tables of results can be downloaded using the following links:

- [**Gene summary results**](./file/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/mageck_gene_summary.xlsx)

- [**sgRNA summary results**](./file/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/mageck_sgrna_summary.xlsx)

For file format description, visit [MAGeCK output](https://sourceforge.net/p/mageck/wiki/output/)

### Results by comparison
Download the results by comparison, including gene-set analysis level:

``` r
x <- foreach(i=names(mageck_comp_list), .combine = c) %do% {
  ilower <- tolower(i) %>% gsub(" ", "_",.) 
  file_xlsx <- file.path('./file', rmd_file, paste0(ilower, '_mageck_results.xlsx'))
  paste0("[", i, "](", file_xlsx,")", sep = "")
}
cat(paste(x, collapse = ' - '), "\n\n")
```

[StealTHY KO Syngeneic CMT167](./file/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/stealthy_ko_syngeneic_cmt167_mageck_results.xlsx) - [StealTHY KO Syngeneic LLC1](./file/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/stealthy_ko_syngeneic_llc1_mageck_results.xlsx) - [StealTHY KO Syngeneic CT26](./file/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/stealthy_ko_syngeneic_ct26_mageck_results.xlsx) 


## Differential abundance analysis

``` r
gene_info <- readxl::read_xlsx('./data/resources/selected_gene_list/gene_list_for_figure_5c_additional_models.xlsx')

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


gene_levels <- rev(gene_info$id)


gene_levels_list <- list()
gene_levels_list[['stealthy_ko_syngeneic_llc1']] <- c(
  "Pcna",
  "Amh", "Amhr2",  "Tgfbr2", "Acvr2b", "Gdf15", "Tgfb3", "Tgfbr1",
  "Mertk", "Igf2r", "Angpt2", "Igf1r", "Egfr", "Axl", "Gas6",
  "Plk1","Trim72",
  "Cd30","Ifnlr1","Ifnar1","Trail","Ifne","Lta","Ltbr","Tnfr2","Light", 
  "Ccl25","Cxcl9","Ccrl2","Cxcr5","Ccr6","Ccr4","Relt", "Il18ra1", "Il13ra2", "Il13ra1", "Il18","Cxcl4","Ccl4","Cxcl14"
  )

gene_levels_list[['stealthy_ko_syngeneic_ct26']] <- c(
  "Pcna",
  'Amhr2', 'Amh', 'Gdf15', 'Tgfbr2', 'Gdf9', 'Acvr1b', 'Tgfbr1',
  'Axl', 'Egfr', 'Erbb3', 'Igf2r', 'Mertk', 'Gas6', 'Met',
  "Plk1","Trim72",
  "Cd30","Ifnlr1","Ifnar1","Trail","Ifne","Lta","Ltbr","Tnfr2","Light", 
  "Ccl25","Cxcl9","Ccrl2","Cxcr5","Ccr6","Ccr4","Relt", "Il18ra1", "Il13ra2", "Il13ra1", "Il18","Cxcl4","Ccl4","Cxcl14"
  )

gene_levels_list[['stealthy_ko_syngeneic_cmt167']]  <- c(
  "Pcna",
  'Amhr2', 'Amh', 'Gdf9', 'Acvr1b', 'Acvr1', 'Tgfbr2', 'Tgfbr1',
  'Met', 'Igf1r', 'Igf2r', 'Egfr', 'Mertk', 'Axl', 'Angpt2',
  "Plk1","Trim72",
  "Cd30","Ifnlr1","Ifnar1","Trail","Ifne","Lta","Ltbr","Tnfr2","Light",
  "Ccl25","Cxcl9","Ccrl2","Cxcr5","Ccr6","Ccr4","Relt", "Il18ra1", "Il13ra2", "Il13ra1", "Il18","Cxcl4","Ccl4","Cxcl14"
  )
```

### Supplementary figure 4M: Differential abundance of selected genes


Dot plots highlighting hits identified following the strategy depicted in Fig. 4A upon transplantation in immunocompetent hosts of metastatic lung carcinoma cells (LLC1, C167) and colon carcinoma cells (CT26) previously subjected to StealTHY KO; representative genes from two recurrently identified categories are shown in comparison across screening methods. The fold-change values indicate differential sgRNA abundance between primary tumors and lung metastases; n=5 each; RTK, receptor tyrosine kinase; Ca, carcinoma; FDR, false discovery rate; Tnf, tumor necrosis factor; Ifn, interferon.


``` r
use_comp <- names(gene_summ)

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

xlfc <- foreach(i=gene_summ, .combine = rbind) %do% {
  i %>%
    ReadRRA(score = 'lfc') %>%
    dplyr::rename(LFC = Score) %>%
    left_join(gene_info, by = c('id' = 'id')) %>%
    filter(!is.na(gene_set)) %>%
    pull(LFC) %>%
    abs %>% max
}

i <- use_comp[1]
for(i in use_comp)  {
  cat("####", mageck_comp_list_names[i], "\n\n")
  
  gene_levels <-  gene_levels_list[[i]]
  
  x <- gene_summ[[i]] %>% 
    ReadRRA(score = 'lfc') %>% 
    dplyr::rename(LFC = Score) %>% 
    left_join(gene_info, by = c('id' = 'id')) %>% 
    filter(!is.na(gene_set)) %>% 
    filter(gene_name %in% gene_levels) %>% 
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
  size_breaks <- c(0, 0.25, 0.5, ceiling(maxFDR))
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
  
  res <- res +
    scale_x_continuous(
      expand = expansion(add = c(1, 1)), 
      limits = xlim
    )
  
  print(res)
  cat('\n\n')
}
```

#### StealTHY KO Syngeneic CMT167 

<img src="figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-1.png" style="display: block; margin: auto;" />

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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-1.png" target="_blank">a6972e9</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

#### StealTHY KO Syngeneic CT26 

<img src="figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-2.png" style="display: block; margin: auto;" />

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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-2.png" target="_blank">a6972e9</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  

#### StealTHY KO Syngeneic LLC1 

<img src="figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-3.png" style="display: block; margin: auto;" />

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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/da-selected-genes-dotplot-3.png" target="_blank">a6972e9</a></td>
  <td>Francesc Castro-Giner</td>
  <td>2025-07-01</td>
  </tr>
  </tbody>
  </table>
  </div>
  </div>
  


## Supplementary figure 4N: Muller plots of selected genes

Muller plots depicting the relative frequency of unique gene KOs in matched primary lung tumors versus synchronous liver metastasis samples from syngeneic recipients transplanted with C167 lung carcinoma cells, previously transduced with the mouse interactome library; n=5 (two representative samples are shown).


``` r
x <- colData(se) %>% data.frame

mageck_comp_list <- list(
  `StealTHY KO Syngeneic CMT167` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'CMT167-KO-Thy1-C57BL6' & sample_type == 'Primary tumor' ) %>% pull(sample_alias),
    `Liver mets` = x %>% filter(vector_mmodel == 'CMT167-KO-Thy1-C57BL6' & sample_type == 'Liver mets' ) %>% pull(sample_alias)
    ),
  `StealTHY KO Syngeneic LLC1` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'LLC1-KO-Thy1-C57BL6' & sample_type == 'Primary tumor' ) %>% pull(sample_alias),
    `Liver mets` = x %>% filter(vector_mmodel == 'LLC1-KO-Thy1-C57BL6' & sample_type == 'Liver mets' ) %>% pull(sample_alias)
    ),
  `StealTHY KO Syngeneic CT26` = list(
    `Primary tumor` = x %>% filter(vector_mmodel == 'CT26-KO-Thy1-BALB' & sample_type == 'Primary tumor' ) %>% pull(sample_alias),
    `Liver mets` = x %>% filter(vector_mmodel == 'CT26-KO-Thy1-BALB' & sample_type == 'Lung mets' ) %>% pull(sample_alias)
    )
)


comp_list_muller <- list(
  `Comparison 1` = list(
    `CMT167 Syngeneic Thy1.1 4` = c('C167_StealTHY_Tum_4',	'C167_StealTHY_Met_4'),
    `CMT167 Syngeneic Thy1.1 5` = c('C167_StealTHY_Tum_5',	'C167_StealTHY_Met_5')
  )
)
genes_list_1 <- gene_info %>% 
  filter(gene_set %in% c("TGF-β family related", "RTK signaling", "Positive controls")) %>% 
  pull(id)
genes_list_2 <- gene_info %>% 
  filter(!gene_set %in% c("TGF-β family related", "RTK signaling", "Positive controls")) %>% 
  pull(id)
```



``` r
cpm_threshold <- 10000
comp_list <- comp_list_muller
use_genes <- c(
  "Egfr",
  'Met',
  'Igf1r',
  'Igf2r',
  'Amhr2',
  'Tgfbr1',
  'Tgfb2',
  'Tgfb3',
  'Gdf9',
  'Mertk',
  'Gas9',
  'Axl',
  'Trim72'
)


# Create DF
use_df <- assay(se, 'cpm') %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('guide') %>% 
  mutate(
    gene = rowData(se)$Gene,
    gene = ifelse(grepl('Non_Target', gene), 'Non_Target', gene)
  ) %>%
  filter(gene %in% use_genes) %>% 
  pivot_longer(-c(guide, gene), names_to = 'sample_alias', values_to = 'cpm')


j <-  names(comp_list)[1]
muller_plots <- foreach(j = names(comp_list)) %do% {
  i <- comp_list[[j]][[1]]
  iname <-  names(comp_list[[j]])[1]
  res <- foreach(iname = names(comp_list[[j]])) %do% {
    i <- comp_list[[j]][[iname]]
    x <- use_df %>% 
      filter(sample_alias %in% i) %>% 
      mutate(
        sample_type = ifelse(sample_alias == i[1], 'Primary tumor', 'Metastasis'),
        sample_type = factor(sample_type, levels = c('Primary tumor', 'Metastasis'))
      )
    
    # Group CPM by gene
    x %<>% 
      group_by(gene, sample_alias, sample_type) %>% 
      summarise(cpm = sum(cpm)) %>% 
      mutate(
        gene = factor(gene, levels = use_genes)
      )
    use_cols <- c(
      'low abundant' = 'grey80',
      colorRampPalette(palette_OkabeIto)(nlevels(x$gene)) %>% set_names(levels(x$gene))
    )
    
    x %>% 
      ggplot( aes(x = sample_type, y = cpm/10000, group = gene, fill = gene)) + 
      geom_area(colour = alpha("white", 0.1), linewidth = 0.08, alpha = 0.8) +
      scale_fill_manual(values = use_cols, guide = guide_legend(ncol = 3)) +
      labs(
        x = '',
        y = expression("CPM x 10"^-4),
        title = iname,
        fill = ''
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      # guides(fill = 'none') +
      theme(
        panel.background = element_rect(fill = "grey90"),
        plot.margin = margin(0, 0.5, 0, 0, "cm"),
        axis.text =  element_text(size=3),
        axis.title =  element_text(size=3),
        plot.title = element_text(size=3, hjust = 0.5),
        legend.text = element_text(size=3)
      )
  }
  names(res) <- names(comp_list[[j]])
  return(res)

}
names(muller_plots) <- names(comp_list)
```



``` r
for (i in names(muller_plots)) {
  cat('###', i, '\n\n')
  wrap_plots(muller_plots[[i]], ncol = 1) %>% print
  cat('\n\n')
}
```

### Comparison 1 

<img src="figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/muller-plots-1.png" style="display: block; margin: auto;" />

  <p>
  <button type="button" class="btn btn-default btn-xs btn-workflowr btn-workflowr-fig"
  data-toggle="collapse" data-target="#fig-muller-plots-1">
  Past versions of muller-plots-1.png
  </button>
  </p>

  <div id="fig-muller-plots-1" class="collapse">
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
  <td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/a6972e9baca3624749a8db52f17092d7d2e8c2e2/docs/figure/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.Rmd/muller-plots-1.png" target="_blank">a6972e9</a></td>
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
 [9] ggridges_0.5.6              openxlsx_4.2.8             
[11] SummarizedExperiment_1.36.0 Biobase_2.66.0             
[13] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[15] IRanges_2.40.1              S4Vectors_0.44.0           
[17] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
[19] matrixStats_1.5.0           diptest_0.77-1             
[21] kableExtra_1.4.0            DT_0.34.0                  
[23] magrittr_2.0.3              foreach_1.5.2              
[25] knitr_1.50                  lubridate_1.9.4            
[27] forcats_1.0.0               stringr_1.5.1              
[29] dplyr_1.1.4                 purrr_1.1.0                
[31] readr_2.1.5                 tidyr_1.3.1                
[33] tibble_3.3.0                ggplot2_3.5.2              
[35] tidyverse_2.0.0             workflowr_1.7.1            

loaded via a namespace (and not attached):
  [1] splines_4.4.3           later_1.4.4             ggplotify_0.1.2        
  [4] bitops_1.0-9            filelock_1.0.3          cellranger_1.1.0       
  [7] R.oo_1.27.1             graph_1.84.1            XML_3.99-0.19          
 [10] lifecycle_1.0.4         httr2_1.2.1             rstatix_0.7.2          
 [13] rprojroot_2.1.0         processx_3.8.6          lattice_0.22-7         
 [16] crosstalk_1.2.2         backports_1.5.0         sass_0.4.10            
 [19] rmarkdown_2.29          jquerylib_0.1.4         yaml_2.3.10            
 [22] ggtangle_0.0.7          httpuv_1.6.16           zip_2.3.3              
 [25] depmap_1.20.0           DBI_1.2.3               RColorBrewer_1.1-3     
 [28] abind_1.4-8             zlibbioc_1.52.0         R.utils_2.13.0         
 [31] RCurl_1.98-1.17         yulab.utils_0.2.0       rappdirs_0.3.3         
 [34] git2r_0.36.2            GenomeInfoDbData_1.2.13 enrichplot_1.26.6      
 [37] ggrepel_0.9.6           tidytree_0.4.6          svglite_2.2.1          
 [40] codetools_0.2-20        DelayedArray_0.32.0     DOSE_4.0.1             
 [43] xml2_1.3.8              tidyselect_1.2.1        aplot_0.2.8            
 [46] UCSC.utils_1.2.0        farver_2.1.2            BiocFileCache_2.14.0   
 [49] pathview_1.46.0         jsonlite_2.0.0          Formula_1.2-5          
 [52] iterators_1.0.14        systemfonts_1.2.3       tools_4.4.3            
 [55] treeio_1.30.0           Rcpp_1.1.0              glue_1.8.0             
 [58] gridExtra_2.3           SparseArray_1.6.2       xfun_0.53              
 [61] qvalue_2.38.0           withr_3.0.2             BiocManager_1.30.26    
 [64] fastmap_1.2.0           callr_3.7.6             digest_0.6.37          
 [67] gridGraphics_0.5-1      timechange_0.3.0        R6_2.6.1               
 [70] textshaping_1.0.1       GO.db_3.20.0            dichromat_2.0-0.1      
 [73] RSQLite_2.4.3           R.methodsS3_1.8.2       generics_0.1.4         
 [76] data.table_1.17.8       httr_1.4.7              htmlwidgets_1.6.4      
 [79] S4Arrays_1.6.0          whisker_0.4.1           pkgconfig_2.0.3        
 [82] gtable_0.3.6            blob_1.2.4              XVector_0.46.0         
 [85] clusterProfiler_4.14.6  htmltools_0.5.8.1       carData_3.0-5          
 [88] fgsea_1.32.4            scales_1.4.0            png_0.1-8              
 [91] ggfun_0.2.0             rstudioapi_0.17.1       tzdb_0.5.0             
 [94] reshape2_1.4.4          nlme_3.1-168            curl_7.0.0             
 [97] org.Hs.eg.db_3.20.0     cachem_1.1.0            BiocVersion_3.20.0     
[100] parallel_4.4.3          vipor_0.4.7             AnnotationDbi_1.68.0   
[103] pillar_1.11.0           grid_4.4.3              vctrs_0.6.5            
[106] promises_1.3.3          car_3.1-3               dbplyr_2.5.0           
[109] beeswarm_0.4.0          Rgraphviz_2.50.0        evaluate_1.0.5         
[112] KEGGgraph_1.66.0        cli_3.6.5               compiler_4.4.3         
[115] rlang_1.1.6             crayon_1.5.3            ggsignif_0.6.4         
[118] labeling_0.4.3          ps_1.9.1                getPass_0.2-4          
[121] plyr_1.8.9              fs_1.6.6                stringi_1.8.7          
[124] viridisLite_0.4.2       BiocParallel_1.40.2     Biostrings_2.74.1      
[127] lazyeval_0.2.2          GOSemSim_2.32.0         Matrix_1.7-3           
[130] ExperimentHub_2.14.0    hms_1.1.3               bit64_4.6.0-1          
[133] KEGGREST_1.46.0         AnnotationHub_3.14.0    igraph_2.1.4           
[136] broom_1.0.9             memoise_2.0.1           bslib_0.9.0            
[139] ggtree_3.14.0           fastmatch_1.1-6         bit_4.6.0              
[142] readxl_1.4.5            gson_0.1.0              ape_5.8-1              
</div>
