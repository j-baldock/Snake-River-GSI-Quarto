library(quarto)
library(knitr)

# Render all
quarto::quarto_render(output_format = "html")

# Render 'chapt_name' chpater only
quarto::quarto_render("index.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Baseline Relatedness/BaselineRelatedness.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Baseline Testing/BaselineDataTesting.qmd", 
                      cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("BaselineFst/BaselineFst.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("GSI Analysis/GSIAnalysis.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Landscape Covariates/Watershed Delineation/WatershedDelineation.qmd", 
                      cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Landscape Covariates/Flowline Distance/FlowlineDist.qmd", 
                      cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Landscape Covariates/Groundwater/GroundwaterIndex.qmd", 
                      cache_refresh = TRUE, # default is FALSE
                      output_format = "html")
