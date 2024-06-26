# Scater palettes.
# Obtained from : https://github.com/alanocallaghan/scater/blob/devel/R/plot_colours.R
get_scater_palette <- function(palette_name)
  # Function to define colour palettes.
{
  switch(palette_name,
         tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                       "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                       "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                       "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
         tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                             "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                             "#CDCC5D", "#6DCCDA"),
         colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                          "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                          "#FFBC79", "#CFCFCF"),
         colourblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
                           "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
                           "#FFBC79", "#CFCFCF"),
         trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
                          "#FFC156", "#69B764", "#F26C64", "#FFDD71",
                          "#9FCD99"),
         purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
                          "#5F5A41", "#B4B19B", "#995688", "#D898BA",
                          "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
         bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
                       "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
                       "#BD0A36", "#F4737A"),
         greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
                           "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
                           "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
         cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
                    "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
                    "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
                    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
         okabeito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

  )
}


# Specific palettes
palette_cell_type <- c(
  monocyte_like = '#F0E442',
  monocyte = '#F0E442',
  Monocyte = '#F0E442',
  neutrophil_like = '#E69F00',
  neutrophil = '#E69F00',
  Neutrophil = '#E69F00',
  nk_cells = '#D55E00',
  `NK cells` = '#D55E00',
  `NK` = '#D55E00',
  t_cells  = '#0072B2',
  `T cells`  = '#0072B2',
  `T cell`  = '#0072B2')

palette_donor <- c(
  `4T1` = '#B5C8E2',
  MVT1 = '#CDCC5D',
  Py2T = '#F4737A'
  )


palette_vector <- c(
  `Thy1`  = '#F7EC13',
  `Puro`  = '#58595B',
  `dCas9+Puro` = 'black'
)

palette_vector_line <- c(
  `Thy1` = 'black',
  `Puro`  = 'white',
  `dCas9+Puro` = 'white'
)


palette_mmodel <- c(
  NSG = '#58595B',
  Syngeneic = '#F7EC13'
)

palette_mmodel_line <- c(
  NSG = 'white',
  Syngeneic = 'black'
)


palette_mmodel_line <- c(
  `StealTHY syngeneic` = 'black',
  `StealTHY NSG` = 'white'
)

palette_vector_mmodel <- c(
  `StealTHY syngeneic`  = '#F7EC13',
  `StealTHY NSG`  = '#58595B'
)


palette_immunogen_burden <- c(
  `HuTHY1`  = '#F7EC13',
  `GFP`  = '#58595B',
  `Cas9-HuTHY1`  = '#58595B',
  `Cas9-GFP` = 'black'
)


palette_hu_stealthy_ko <- c(
  `huNSG`  = '#F7EC13',
  `NSG`  = '#58595B'
)



