![alt text](./Shiny/www/logo_intervista_new.png "Inter-ViSTA (portable)")

Inter-ViSTA (**Inter**action **Vi**sualization in **S**pace and **T**ime **A**nalysis) is a computational platform designed by the Cristea Lab at Princeton University to integrate spatial and temporal proteome and protein interaction information from user derived experiments with the wealth of knowledge in existing databases. 

Inter-ViSTA is available as both a web-based platform ([Inter-ViSTA online](http://intervista.princeton.edu)), as a portable program when downloaded from the [Inter-ViSTA github page](https://github.com/cristealab/Inter-ViSTA), or as a [Docker image](https://hub.docker.com/repository/docker/cristealab/intervista).

## Why use Inter-ViSTA?
Inter-ViSTA utilizes cutting-edge graph visualization algorithms and automatic programmatic access to protein databases to deliver an intuitive and user-friendly data visualization platform. The platform enables users to discover underlying patterns in multiple interaction networks across conditions using an interactive interface with dynamic network visuals and exploratory quantitative analysis. Users can integrate interaction networks of one or many baits, include subcellular localization information and functional annotations to modify networks, analyze dynamic quantitative properties of the network, and identify properties of proteins shared between multiple baits - all with their own data, with few manual steps.

## Which branch of the Inter-ViSTA GitHub repository should I use?
We currently maintain 4 slightly different branches of the Inter-ViSTA application on GitHub:

1. Base (general use) - This branch contains only the core Inter-ViSTA application code. For it to function properly, you must have both R and Python (Anaconda) installed on your local machine. See *What if I already have R and Anaconda on my machine?* below for more information on specific dependencies.
2. Portable (general use) - In addition to the core code, this branch contains a portable version of RStudio and the necessary R libraries to run Inter-ViSTA. However, it still requires a local Python installation and a properly-configured conda environment to run properly. See *System requirements and dependencies* and *What if I don't have Python on my computer?*
3. Docker (for development/deployment only) - Core code + necessary files to generate the Docker image found at https://hub.docker.com/repository/docker/cristealab/intervista
4. shinyapps.io (for development/deployment only) - Core code modified to run on the system configuration required by shinyapps.io

## Running an Inter-ViSTA analysis
Inter-ViSTA is designed to be compatible with a wide variety of interaction-based proteomics approaches (AP-MS, APEX, BioID, etc.). See the "Instructions" tab of the application to learn how to format and upload your own data for visualization and analysis with Inter-ViSTA.

Additionally, we have provided an instructional video illustrating the capabilities of Inter-ViSTA below:

[![Alt text](https://img.youtube.com/vi/9d80gxNqKTw/0.jpg)](https://www.youtube.com/watch?v=9d80gxNqKTw)

## System requirements and dependencies
*Note: The following dependencies apply only to the portable version of Inter-ViSTA. [Inter-ViSTA online](http://intervista.princeton.edu) only requires access to a web browser and the internet.*

Inter-ViSTA combines the stunning network vizualization capacities of R with the API querying abilities of Python and therefore requires both for full functionality. Conveniently, we have packaged Inter-ViSTA with a portable version of R and an accompanying web browser so that the application can be run locally with as few software dependencies as possible. Therefore, the only additional software required to run Inter-ViSTA is Anaconda Python 3.7 and several underlying Python packages.

### What if I don't have Python on my computer?
If you do not have Python installed on your system, we recommend downloading it via [Anaconda](https://www.anaconda.com). After installation, run the following command within an Anaconda prompt to set up a conda environment with the packages required to run Inter-ViSTA.
```
conda create -n intervista biopython requests pandas
```
Then simply download and extract the Inter-ViSTA.zip folder from the main github page and double-click "run.vbs" to run the application.

### What if I already have R and Anaconda on my machine?
If R and Anaconda are already installed, you will only need the files in Inter-ViSTA's ["Shiny" folder](https://github.com/cristealab/Inter-ViSTA/tree/private/Shiny) to run the application.

Prior to initializing the application, ensure that the R dependencies are fulfilled by running the following in the R console:
```
install.packages(c("DT",
                   "shiny",
                   "htmlwidgets",
                   "shinyBS",
                   "shinyjs",
                   "zip",
                   "ggplot2",
                   "gridExtra",
                   "reticulate",
                   "networkDynamic",
                   "RColorBrewer",
                   "ndtv",
                   "tsna"
                   ))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("topGO")
```
 Additonally, ensure that Biopython, requests, and pandas are installed on your base Anaconda environment or in a conda environment named "intervista".
 
 ## Citing Inter-ViSTA
 Inter-ViSTA is provided under the MIT license. Specific publication details will be updated in the future.
