---
title: "StealTHY eludes CRISPR immunogenicity to expose concealed metastatic regulators"
author: "Francesc Castro-Giner"
site: workflowr::wflow_site
output:
  workflowr::wflow_html:
    toc: true
editor_options:
  chunk_output_type: console
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
  2
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
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongRepositoryversionstrongahrefhttpsgithubcomTheAcetoLabsainistealTHYtree54b3c4df86314709e556aa3999fc0bbf4226f47dtargetblank54b3c4da">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Repository version:</strong> <a href="https://github.com/TheAcetoLab/saini-stealTHY/tree/54b3c4df86314709e556aa3999fc0bbf4226f47d" target="_blank">54b3c4d</a>
</a>
</p>
</div>
<div id="strongRepositoryversionstrongahrefhttpsgithubcomTheAcetoLabsainistealTHYtree54b3c4df86314709e556aa3999fc0bbf4226f47dtargetblank54b3c4da" class="panel-collapse collapse">
<div class="panel-body">
  
<p>
Great! You are using Git for version control. Tracking code development and
connecting the code version to the results is critical for reproducibility.
</p>

<p>
The results in this page were generated with repository version <a href="https://github.com/TheAcetoLab/saini-stealTHY/tree/54b3c4df86314709e556aa3999fc0bbf4226f47d" target="_blank">54b3c4d</a>.
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
	Untracked:  analysis/about.md
	Untracked:  analysis/crispr-hsapiens_2180_sgRNA.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-StealTHY.md
	Untracked:  analysis/crispr-mm_2215_sgRNA-clonality.md
	Untracked:  analysis/index.md
	Untracked:  analysis/license.md
	Untracked:  analysis/rnaseq-amhr2_ko-deg.md
	Untracked:  analysis/rnaseq-tumor-bulk.md
	Untracked:  analysis/rnaseq-tumor-facs-cancer.md
	Untracked:  analysis/rnaseq-tumor-facs-immune.md
	Untracked:  analysis/style.css
	Untracked:  analysis/tcga_immune_infiltrate_amhr2_ko_signatures.md
	Untracked:  analysis/tcga_survival_amhr2_ko_signatures.md
	Untracked:  analysis/tcga_survival_crispr_hits_genes.md
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
	Deleted:    update_workflowr.R

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
to the R Markdown (<code>analysis/index.Rmd</code>) and HTML (<code>docs/index.html</code>)
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
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/54b3c4df86314709e556aa3999fc0bbf4226f47d/analysis/index.Rmd" target="_blank">54b3c4d</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Update index, license and about</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/5e831e20de8357cadc402a20a98aa789b537a814/docs/index.html" target="_blank">5e831e2</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/c44a44f06465066cb175aacc9cfe2cbe6376900a/analysis/index.Rmd" target="_blank">c44a44f</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Update index, license and about</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/595c4b0d5d16f835cd86f5fe463b47ce59cd3cea/docs/index.html" target="_blank">595c4b0</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/aa04e35a434c6815631b0107424f4798d3fd697b/analysis/index.Rmd" target="_blank">aa04e35</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Update index, license and about</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/4921153ca8c84dd525d308e4183929b6d9909a0f/docs/index.html" target="_blank">4921153</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/455fa9121fcca51363cf639f7ee170fb108e34c1/analysis/index.Rmd" target="_blank">455fa91</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Prepare files for publication</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/13b89b011934eb3bd6ea30828c0fe8a0c0cd5a92/docs/index.html" target="_blank">13b89b0</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-17</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/5c3fb91479e6eb57971b4b2e04b5d7c95d632e4a/analysis/index.Rmd" target="_blank">5c3fb91</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-17</td>
<td>Start workflowr project.</td>
</tr>
</tbody>
</table>
</div>

<hr>
</div>
</div>
</div>




## Publication

Saini, M., Castro-Giner, F., Hotz, A., Sznurkowska, M.K., Nüesch, M., Cuenot, F.M., Budinjas, S., Bilfeld, G., Ildız, E.S., Strittmatter, K., Krol, I., Diamantopoulou, Z., Paasinen-Sohns, A., Waldmeier, M., Cássio, R., Kreutzer, S., Kontarakis, Z., Gvozdenovic, A., Aceto, N. StealTHY eludes CRISPR immunogenicity to expose concealed metastatic regulators (2025)


## Abstract

CRISPR screens have become standard gene discovery platforms in various contexts, including cancer. Yet, commonly available CRISPR-Cas9 tools are increasingly recognized as unfit for in vivo investigations in immunocompetent contexts, due to broad immunogenicity of bacterial nucleases and reporters. Here, we show how conventional CRISPR screens in tumor grafts are systematically jeopardized by immunoediting in syngeneic and humanized immunocompetent hosts, resulting in iatrogenic clonal dropouts and ultimately compromising target identification. To resolve this, we present StealTHY, an immunogen-free CRISPR platform compatible with virtually all immunocompetent designs, enabling preservation of clonal architecture and exposing previously concealed cancer vulnerabilities. Among these, we identify the AMH-AMHR2 axis as a formerly unappreciated metastasis target. Thus, with StealTHY, we provide a new resource to expand the applicability of CRISPR screens to immunocompetent models including humanized tumor grafts, revealing metastasis regulators of therapeutic relevance. 

## Data pre-processing

Raw data is available at Gene Expression Omnibus (GEO, NCBI), under the accession codes [GSE268044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268044), [GSE268045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268045) and [GSE268588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268588). Data processing is computationally expensive and is not covered in this repository. We provide description of the data pre-processing workflow together with software version in the original publication. Processed data, large result files, additional functions, references and metadata are were archived at [![DOI](https://zenodo.org/doi/10.5281/zenodo.11208099.svg)](https://doi.org/10.5281/zenodo.11208099)

##  Data and code availability

To reproduce our analysis, first clone source code from the [GitHub repository](https://github.com/TheAcetoLab/saini-stealTHY). This repository is also archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11288553.svg)](https://doi.org/10.5281/11288553)

    git clone https://github.com/TheAcetoLab/saini-stealTHY

Next, download processed data and output data deposited in [Zenodo](https://doi.org/10.5281/zenodo.11374563) into the cloned project folder and untar the files.

    for file in *.tar.gz; do tar xzvf "${file}" && rm "${file}"; done
    

## Reproducibility

The results form our analyses are listed below in webpage format. They were generated from R Markdown documents deposited in the [GitHub repository](https://github.com/TheAcetoLab/saini-stealTHY). The workflow of the analysis was created using the [workflowr](https://cran.r-project.org/web/packages/workflowr/index.html) R package and can be reproduced in its totality using [workflowr](https://cran.r-project.org/web/packages/workflowr/index.html) [wflow_build](https://jdblischak.github.io/workflowrBeta/reference/wflow_build.html) command after the installation of the proper R packages. Session info, including R and package versions, was automatically included at the end of each analysis file.

Files containing pre-computed results from differential expression or gene-set enrichment analyses were deposited in [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11374563.svg)](https://doi.org/10.5281/zenodo.11374563). In order to generate those files again change the option `eval = FALSE` to `eval = TRUE` in the specific code chunk from the R Markdown file.

## Analyses
-   [In-vivo CRISPR screening in mouse: sgRNA distribution](crispr-mm_2215_sgRNA-clonality.html)
-   [RNA-seq of primary tumor transcriptomes from syngeneic hosts based on Cas9 delivery TrCD8a status](rnaseq-tumor-bulk.html)
-   [Restoration of metastatic gene discovery by StealTHY CRISPR screening](crispr-mm_2215_sgRNA-StealTHY.html)
-   [CRISPR screening transplantation in immunocompetent hosts of metastatic lung carcinoma cells (LLC1, C167) and colon carcinoma cells (CT26) previously subjected to StealTHY KO](crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.html)
-   [CRISPR screening in humanized mice by StealTHY](crispr-hsapiens_2180_sgRNA.html)
-   [Survival analysis of gene-signatures identified in the CRISPR screening in humanized mice by StealTHY](tcga_survival_crispr_hits_signatures.html)
-   [RNA-sequencing of primary tumors from syngeneic FVB/N orthotopically transplanted with Py2T carrying StealTHY KO of Amhr2 or a non-targeting sgRNA](rnaseq-amhr2_ko-deg.html)
-   [RNA-seq of FACs sorted CD3+ TIL cells comparing Cas9-Thy1(A) versus Thy1(A) transduction](rnaseq-tumor-facs-immune.html)
-   [RNA-seq of FACs sorted carcinoma cells comparing Cas9-Thy1(A) versus Thy1(A) transduction](rnaseq-tumor-facs-cancer.html)
-   [Correlation of *AMH* signatures with immune cell fractions in TCGA](tcga_immune_infiltrate_amhr2_ko_signatures.html)
-   [Survival analysis of *AMH* signatures in TCGA cohorts](tcga_survival_amhr2_ko_signatures.html)
-   [Survival analysis of *AMH*, *AMHR2* and *GDF9* expression in TCGA cohorts](tcga_survival_crispr_hits_genes.html)



## Main figures

### Figure 2
-   [Figure 2B](crispr-mm_2215_sgRNA-clonality.html#31_Figure_2B:_sgRNA_distribution_in_mice_transplanted_with_4T1_cells)
-   [Figure 2C](crispr-mm_2215_sgRNA-clonality.html#5_Figure_2C:_Representative_Müller_plots_of_4T1_in_syngeneic_recipients)
    
### Figure 3
-   [Figure 3D](rnaseq-tumor-bulk.html#Figure_3D_:_qRT-PCR_in_TrCD8a+_and_TrCD8a-)
-   [Figure 3H](rnaseq-tumor-bulk.html#Figure_3H:_Sample_distances_in_MVT1_models)

### Figure 4
-   [Figure 4D](crispr-mm_2215_sgRNA-StealTHY.html#sgRNA_distribution2)
-   [Figure 4C](crispr-mm_2215_sgRNA-StealTHY.html#Figure_4C_and_Supplementary_Figure_4F)
-   [Figure 4D](crispr-mm_2215_sgRNA-StealTHY.html#Figure_4D:_Proportion_of_non-expressed_hits)

### Figure 5   
-   [Figure 5H](crispr-hsapiens_2180_sgRNA.html#31_Figure_5H:_Muller_plots_for_huNSGs)

### Figure 6    
-   [Figure 6B](crispr-hsapiens_2180_sgRNA.html#sgRNA_distribution_by_condition_Primary_tumor)
-   [Figure 6C](crispr-hsapiens_2180_sgRNA.html#Figure_6C:_Volcano_plots_showing_genes_Lungs_over_Primary_Tumor_in_StealTHY_KO)
-   [Figure 6D](crispr-hsapiens_2180_sgRNA.html#Figure_6D:_Distribution_of_the_fold-change_of_Lungs_over_Primary_Tumor_in_huNSG_and_NSG_host)
-   [Figure 6E](tcga_survival_crispr_hits_signatures.html#Figure_6F:_Forest_plot_by_study)
-   [Figure 6F](tcga_survival_crispr_hits_signatures.html#Figure_6E:_Kaplan-Meier_plots)

### Figure 7    
-   [Figure 7F](rnaseq-amhr2_ko-deg.html#__531_Comparison_InVivo_AMHR2_KO_+_AMHR2_KO_culture–over–InVivo_Control_sgRNA_+_AMHR2_Control_culture_)
-   [Figure 7G](tcga_immune_infiltrate_amhr2_ko_signatures.html#4_Figure_7G:_Correlation_of_AMH_signatures_with_immune_cell_fractions_across_TCGA_cohorts_(Heatmaps))
-   [Figure 7H: Left panel](tcga_survival_amhr2_ko_signatures.html#Figure_7H_Right:_Pooled_Kaplan-Meier_plots_-_Optimal_cut-off)
-   [Figure 7H: Right panel](tcga_survival_amhr2_ko_signatures.html#Figure_7H_Left:_Forest_plot_by_study)
      

    
## Supplementary figures

### Supplementary Figure 2
-   [Supplementary Figure 2B](crispr-mm_2215_sgRNA-clonality.html#32_Supplementary_Figure_2B:_sgRNA_distribution_in_mice_transplanted_with_Py2T_and_MVT1_cells)
-   [Supplementary Figure 2C-left](crispr-mm_2215_sgRNA-clonality.html#41_Supplementary_Figure_2C_left_:_K-S_test_table)
-   [Supplementary Figure 2C-right](crispr-mm_2215_sgRNA-clonality.html#42_Supplementary_Figure_2C_right_:_Empirical_Cumulative_Density_Function_(ECDF)_plots)
-   [Supplementary Figure 2D: Lung cancer in dCas9 + Puro (left-upper panel)](crispr-mm_2215_sgRNA-clonality.html#21_Supplementary_Figure_2D_(left_upper))
-   [Supplementary Figure 2D: Lung cancer in Thy1 (right-upper panel)](crispr-mm_2215_sgRNA-clonality.html#22_Supplementary_Figure_2D_(right_upper))
-   [Supplementary Figure 2D: Mammary cancer in Puro (left-lower panel)](crispr-mm_2215_sgRNA-clonality.html#Supplementary_Figure_2D:_Metastasis_-_Puro_(left_lower))
-   [Supplementary Figure 2D: Mammary cancer in Thy1 (right-lower panel)](crispr-mm_2215_sgRNA-clonality.html#Supplementary_Figure_2D:_Metastasis_-_Thy1_(right_lower))
-   [Supplementary Figure 2F: Right panel](crispr-mm_2215_sgRNA-clonality.html#23_Supplementary_Figure_2F)
-   [Supplementary Figure 2G: Upper panel](rnaseq-tumor-bulk.html#Supplementary_Figure_2G_top:_Differential_expressed_genes_in_dCas9-Thy1_over_Thy1_transduction)
-   [Supplementary Figure 2G: Lower panel](rnaseq-tumor-bulk.html#Supplementary_Figure_2G_bottom:_Representative_gene-sets_enriched_in_dCas9-Thy1_over_Thy1_transduction)
-   [Supplementary Figure 2H: CD3+ TIL upper panel](rnaseq-tumor-facs-immune.html#41_Supplementary_Figure_2H_top_left:_Differential_expressed_genes_in_dCas9-Thy1_over_Thy1_transduction_in_sorted_CD3-positive_TILs)
-   [Supplementary Figure 2H: CD3+ TIL lower panel](rnaseq-tumor-facs-immune.html#51_Supplementary_Figure_2H_bottom_left:_Representative_gene-sets_enriched_in_dCas9-Thy1_over_Thy1_transduction) 
-   [Supplementary Figure 2H: Carcinoma cells upper panel](rnaseq-tumor-facs-cancer.html#41_Supplementary_Figure_2H_top_right:_Differential_expressed_genes_in_dCas9-Thy1_over_Thy1_transduction_in_sorted_carcinoma_cells)
-   [Supplementary Figure 2H: Carcinoma cells lower panel](rnaseq-tumor-facs-cancer.html#51_Supplementary_Figure_2H_bottom_right:_Representative_gene-sets_enriched_in_dCas9-Thy1_over_Thy1_transduction)  

### Supplementary Figure 3   
-   [Supplementary Figure 3E: Py2T1 (left panel)](rnaseq-tumor-bulk.html#Supplementary_figure_3E_left:_Sample_distances_in_Py2T_model)
-   [Supplementary Figure 3E: 4T1 (right panel)](rnaseq-tumor-bulk.html#Supplementary_figure_3E_right:_Sample_distances_in_4T1_model)

### Supplementary Figure 4  
-   [Supplementary Figure 4B](crispr-mm_2215_sgRNA-StealTHY.html#Supplementary_figure_4B:_sgRNA_distribution_in_mice_transplanted_with_4T1_cells_previously_subjected_to_StealTHY_KO)
-   [Supplementary Figure 4C](crispr-mm_2215_sgRNA-clonality.html#7_Supplementary_Figure_4C)
-   [Supplementary Figure 4D](crispr-mm_2215_sgRNA-StealTHY.html#Supplementary_Figure_4D_:_K-S_test_table)
-   [Supplementary Figure 4E](crispr-mm_2215_sgRNA-StealTHY.html#Supplementary_Figure_4E)
-   [Supplementary Figure 4F](crispr-mm_2215_sgRNA-StealTHY.html#Figure_4C_and_Supplementary_Figure_4F)
-   [Supplementary Figure 4G](crispr-mm_2215_sgRNA-StealTHY.html#Supplementary_figure_4G:_Proportion_of_non-expressed_hits)
-   [Supplementary Figure 4L](crispr-mm_2215_sgRNA-clonality.html#24_Supplementary_Figure_4L)
-   [Supplementary Figure 4M](crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.html#Supplementary_figure_4M:_Differential_abundance_of_selected_genes)
-   [Supplementary Figure 4N](crispr-mm_2215_sgRNA-CMT167-LLC1-CT26.html#Supplementary_figure_4N:_Muller_plots_of_selected_genes)
    
### Supplementary Figure 6
-   [Supplementary Figure 6D](crispr-hsapiens_2180_sgRNA.html#Supplementary_Figure_6D:_Fold-change_of_Lungs_over_Primary_Tumor_in_huNSG_and_NSG_host)
-   [Supplementary Figure 6L](tcga_survival_crispr_hits_genes.html#Supplementary_Figure_6L:_Forest_plot_by_study)
-   [Supplementary Figure 6E](tcga_survival_crispr_hits_genes.html#Supplementary_Figure_6E:_Kaplan-Meier_plots)

### Supplementary Figure 7     
-   [Supplementary Figure 7L](rnaseq-amhr2_ko-deg.html#Supplementary_Figure_7L)
-   [Supplementary Figure 7M: Positive Galunisertib](rnaseq-amhr2_ko-deg.html#Supplementary_Figure_7M_(Galunisertib_positive))
-   [Supplementary Figure 7M: Negative Galunisertib](rnaseq-amhr2_ko-deg.html#Mammary_epithelium_proliferation30)
-   [Supplementary Figure 7P: Epithelial genes](rnaseq-amhr2_ko-deg.html#EMT_regulators_and_targets)
-   [Supplementary Figure 7P: Proliferative genes](rnaseq-amhr2_ko-deg.html#Mammary_epithelium_proliferation)
-   [Supplementary Figure 7U: Negative Galunisertib](rnaseq-amhr2_ko-deg.html#Supplementary_Figure_7U_left)
-   [Supplementary Figure 7U: Positive Galunisertib](rnaseq-amhr2_ko-deg.html#Supplementary_Figure_7U_right_(Galunisertib_positive))



