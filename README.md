# The Associations Between NMR metabolites and BP

This is a generic template, forked from COMHIS project. Clone and rewrite to start a project with the template.

The aim here is to create a model that enables somewhat painless internal reproduction of a project and to ease communication about a project's structure.

**Most of the features here are recommendations, and can obviously be varied on as needed.**

## Repository name

A proposal for unified naming scheme for publication related repositories is as follows: document type_date_project name. An example would be: `article_2019_hume_history_text_reuse`. Date should follow the format YYYYMMDD, with month and day optional, and would probably refer to the projected or actual end date of the project. The whole date -element can also be optional, to be included only if relevant, eg. `article_hume_history_text_reuse` would be equally valid.

## Practices

* **Project overview documentation:**
    * Should reside in the project root in a `README.md` (this file).
    * Should list people involved and their roles in the project.
* **Naming files and folders:**
    * Use all lowercase (except for established standards such as README.md and the .R filename extension).
    * Separate words in file and directory names by underscore: `_`. eg. `south_sea_bubble.R` instead of `south-sea-bubble.R` or `SouthSeaBubble.R`.
* **Structure:**
    * Follow the directory structure laid out below.
    * Include `README.md` in each directory documenting the contents of that directory.
        * This is especially important in data and final code directories.
    * If feasible, to avoid confusion only use single `.gitattributes` and single `.gitignore` file residing in the project root.

## Directory structure

The project repository structured is variation of formats laid out in a few data science project organization articles (see the end of this README). `code` and `output` -directories include `work/` and `final/` -subdirectories. The `work/` -subdirectory is optional, but helps to keep development material separate from the polished and clean end products that should reside in the `final/` directory.

````
rrnmr/
├── README.md               # project overview
├── LICENCE                 # project licence
├── rrnmr-import.Rmd        # R markdown for the data preprocessing
├── articles-importer.R     # Minor R functions for data preprocessing
├──  articles-definitions.R # R functions for variable definitions in data preprocessing
├── rrnmr.Rmd               # R markdown for the analysis (main file)
├── articles-functions.R    # Minor R functions for analysis
├── articles-models.R       # Minor R functions for modeling
├── articles-plots.R        # Minor R functions for plotting 
├── articles-tables.R       # Minor R function for creating tables
````

