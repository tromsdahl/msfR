# MSF R Script Workflows

## PURPOSE
The purpose of this repository is to house the data analysis workflows used by MSF for proteomics, lipidomics, and metabolomics. Using git and github, the workflows will be kept up-to-date and modified with a record of those changes.

## SCOPE
The workflows are separated between **LARGE MOLECULES** and **SMALL MOLECULES**. The bullet lists below detail the workflows currently used.

### LARGE MOLECULES

#### Data Analysis Workflows
- Full Proteomics: Statistical analysis of DIA proteomics from Fragpipe using DIANN
- Phosphoproteomics: Same as above but using the full proteomics to first normalize the phosphopeptides to account for difference in total abundances
  - **NOTE:** This currently doesn't take into account the localization of the phosphorylation.

### SMALL MOLECULES

#### Data Quality Checks
- Data Quality Check (Metabolomics): An R Markdown document that is used to QC LCMS data following peak picking. It should be used prior to the full data analysis.
- Data Quality Check (Lipidomics): The same as above but for lipidomics so that species are separated out by lipid class.

#### Data Analysis Workflows
- Metabolomics: Statistical analysis of metabolomic data following peak picking.
- Lipidomics: Same as above but account for lipid classes

## HOW TO USE

### LARGE MOLECULES
1. Data is first analyzed using FragPipe with DIANN
2. The annotation file and the pg.matrix files are required
3. Follow the individual instructions for the specific workflow being used

### SMALL MOLECULES
1. LCMS data is first integrated to obtain peak areas, height, etc.
2. The data will need to be converted to a .txt file with the metadata added (See wiki for more information)
3. Prior to doing statistical analysis, first perform a Data Quality Check (DQC) to ensure the data was acquired correctly.
4. Once the DQC has passed, then perform the data analysis workflow and follow the individual instructions for that workflow
