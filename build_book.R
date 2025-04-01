library(quarto)
library(knitr)

# Render all
# Delete cache/files for:
#  - SpaceTimeZOID.qmd
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

quarto::quarto_render("Landscape Covariates/Barriers/Barriers.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Landscape Covariates/Landcover/Landcover.qmd", 
                      cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Obj2 Spatiotemporal Var/SpaceTimeZOID.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")

quarto::quarto_render("Obj2 Model Contribution/ModelContributionZOIB.qmd", 
                      #cache_refresh = TRUE, # default is FALSE
                      output_format = "html")
