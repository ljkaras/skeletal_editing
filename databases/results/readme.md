# Data Analysis

### The `results` directory is for storing the results.txt files from the analysis of several databases.

1. The results from each database should be given a subdirectory in the above **`databases`** directory that contains the following:
    * `results.txt`
    * `symmetric_molecules.txt`
2. After a new database is added, `databases.txt` should be updated to include the name of the new subdirectory.
3. After `results` files have been generated for each database, running `final_analysis.py` will generate several dataframes.
    * A dataframe of the number of starting molecules for each heterocycle type will be generated in `products_vs_sms`.
    * A single column version of the above dataframe will be generated in `count_sms_dfs`.
    * Unnormalized dataframes of %new and %common values will be generated in `percent_new_dfs` and `percent_common_dfs`, respecitvely.
    * Unnormalized dataframes of  %new and %common values *relative to the number of SMs* will be generated in `percent_rel_new_dfs` and `percent_rel_common_dfs`, respecitvely.
4. At this stage, the user should navigate to the `count_sms_dfs` and run the `max_common_percentage` script.
    * This will generate several dataframes in `max_common_dfs` for normalization.
5. The user should then navigate to `normalized_common_dfs` and run the `normalize_common_dfs.py` script to generate several **normalized** dataframes, which give unbiased metrics of which transformations generate the fewest known compounds.
6. Running `heatmap_from_csv.py` in any directory will generate heatmaps from each of the dataframes in that folder. 




