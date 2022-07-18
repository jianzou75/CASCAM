<img src="inst/shiny_app/CASCAM_GUI/www/hex-CASCAM.png" width="150">

## CASCAM (Congruence Analysis and Selector of CAncer Models)

CASCAM is an R package for analyzing and selecting the appropriate cancer models (cell lines, PDOs or PDXs) 
based on the patients' tumor information. The figure below shows the framework of our method. 
Several statistics are defined to measure the distance between tumor samples and the cancer models.

<p align="center">
  <img src="inst/shiny_app/CASCAM_GUI/www/Congruence.svg" width="700">
</p>

With the help of CASCAM, users can:
* Identify the genome-wide appropriate cancer models visually and statistically.
* Evaluate the genome-wide pre-selected cancer models in different pathways.
* Explore the similarity between the tumor samples and the cancer model in the pathway topology level.
* Use the above mentioned tools through an R Shiny app interactively.

