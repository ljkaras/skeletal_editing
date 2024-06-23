# Data Analysis

### The `results` directory is for storing the results.txt files from the analysis of several databases.

1. The results from each database should be given a subdirectory in the above **`databases`** directory that contains the following:
    * `results.txt`
    * `symmetric_molecules.txt`
2. After a new database is added, `databases.txt` should be updated to include the name of the new subdirectory.
3. Running `final_analysis.py` will generate several heatmaps from the `results.txt` files for each database. 
    * These plots will be stored in designated subdirectories within the results directory
    * Additionally, the dataframe used to generate each plot will be exported to the appropriate directory within results

