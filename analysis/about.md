---
title: "About"
output:
  workflowr::wflow_html:
    toc: false
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
to the R Markdown (<code>analysis/about.Rmd</code>) and HTML (<code>docs/about.html</code>)
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
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/5e831e20de8357cadc402a20a98aa789b537a814/docs/about.html" target="_blank">5e831e2</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/595c4b0d5d16f835cd86f5fe463b47ce59cd3cea/docs/about.html" target="_blank">595c4b0</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/e544ed7063a82d700e2039eb89e7fbc3cc5b35f9/docs/about.html" target="_blank">e544ed7</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/7bb3fec1f59a54caf60f6f727bd2b08ec5a5dec9/analysis/about.Rmd" target="_blank">7bb3fec</a></td>
<td>Francesc Castro-Giner</td>
<td>2025-10-31</td>
<td>Update license and about</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/TheAcetoLab/saini-stealTHY/13b89b011934eb3bd6ea30828c0fe8a0c0cd5a92/docs/about.html" target="_blank">13b89b0</a></td>
<td>Francesc Castro-Giner</td>
<td>2024-05-17</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/TheAcetoLab/saini-stealTHY/blob/5c3fb91479e6eb57971b4b2e04b5d7c95d632e4a/analysis/about.Rmd" target="_blank">5c3fb91</a></td>
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




## StealTHY eludes CRISPR immunogenicity to expose concealed metastatic regulators

**Abstract**

CRISPR screens have become standard gene discovery platforms in various contexts, including cancer. Yet, commonly available CRISPR-Cas9 tools are increasingly recognized as unfit for in vivo investigations in immunocompetent contexts, due to broad immunogenicity of bacterial nucleases and reporters. Here, we show how conventional CRISPR screens in tumor grafts are systematically jeopardized by immunoediting in syngeneic and humanized immunocompetent hosts, resulting in iatrogenic clonal dropouts and ultimately compromising target identification. To resolve this, we present StealTHY, an immunogen-free CRISPR platform compatible with virtually all immunocompetent designs, enabling preservation of clonal architecture and exposing previously concealed cancer vulnerabilities. Among these, we identify the AMH-AMHR2 axis as a formerly unappreciated metastasis target. Thus, with StealTHY, we provide a new resource to expand the applicability of CRISPR screens to immunocompetent models including humanized tumor grafts, revealing metastasis regulators of therapeutic relevance. 

**Authors**

Massimo Saini<sup>1</sup>, Francesc Castro-Giner<sup>1</sup>, Adriana Hotz<sup>2</sup>, Magdalena K. Sznurkowska<sup>1</sup>, Manuel Nüesch<sup>1</sup>, François M. Cuenot<sup>1</sup>, Selina Budinjas<sup>1</sup> Gilles Bilfeld<sup>1</sup>, Ece Su Ildız<sup>1</sup>, Karin Strittmatter<sup>1</sup>, Ilona Krol<sup>1</sup>, Zoi Diamantopoulou<sup>1</sup>,<sup>3</sup>, Aino Paasinen-Sohns<sup>4</sup>, Maria Waldmeier<sup>1</sup>, Rafaela Cássio<sup>2</sup>, Susanne Kreutzer<sup>2</sup>, Zacharias Kontarakis<sup>2,5</sup>, Ana Gvozdenovic<sup>1</sup>, & Nicola Aceto<sup>1,6,*</sup>

<sup>1</sup> Department of Biology, Institute of Molecular Health Sciences, Swiss Federal Institute of Technology (ETH) Zurich, 8093 Zurich, Switzerland 
<sup>2</sup> Functional Genomics Center Zurich, Swiss Federal Institute of Technology (ETH) Zurich, 8093 Zurich, Switzerland.
<sup>3</sup> Current address: Cancer Research UK Scotland Institute, Garscube Estate, Switchback Road Glasgow, G61 1BD, UK and School of Cancer Sciences, University of Glasgow, Garscube Estate, Switchback Road, Glasgow, G61 1QH, UK
<sup>4</sup> Department of Biomedicine, University of Basel and University Hospital Basel, 4058 Basel, Switzerland
<sup>5</sup> Genome Engineering and Measurement Laboratory (GEML), Swiss Federal Institute of Technology (ETH) Zurich, 8093 Zurich, Switzerland   
<sup>6</sup> Lead Contact

<sup>*</sup>Corresponding author

## Citation


