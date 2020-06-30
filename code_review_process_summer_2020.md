### Code cleaning notes and process MRB beginning 2020-06-30

#### 01_data_cleaning.R
**Doing**:

**To do**:
 + Can the commented out code under `# Make wide version of periphyton using proportion data` be removed?
 + Could consider making all cleaned dataset names lowercase (some are, not all)

**Done**:
 + Changed `filter(Sample_ID ==...)` chain to a single `filter(Sample_ID %in% ...)`
 + Changed `filter(Station == ...)` chain to a single `filter(Station %in% ...)`
 + Linting

#### 02_sewage_indicator_analysis.R
**Doing**:

**To do**:
 + The comment at the top of the script is incomplete

**Done**:
 + Edited some formatting (spacing)
 + Changed to include `"stringsAsFactors" = FALSE` in the `read.csv()` statements. Wasn't a huge issue, but I get nervous when joins start leaving messages that factor levels aren't matching and are being automatically reformatted into strings.
 + Linting

#### 03_community_composition_analysis.R
**Doing**:

**To do**:

**Done**:
 + Changed comment at the top of the script to have a single hashtag
 + Changed to include `"stringsAsFactors" = FALSE` in the `read.csv()` statements, mostly for consistency with other scripts.
 + Some formatting changes
 + Added `ggpubr` to library calls for the `ggarrange()` at the end of the script
 + Removed commented out code from `periphyton_IDW_pop_group_plot` call
 + Linting


 #### 04_fatty_acid_analysis.R
 **Doing**:

 **To do**:
  + Does the "Drapa" replacement in the EFA data_scores have an accidental space before it? (" Drapa spp.")  

 **Done**:
  + Some formatting changes
  + Changed to include `"stringsAsFactors" = FALSE` in the `read.csv()` statements, mostly for consistency with other scripts.
  + Specified `by = c("Site")` for the ppcp_meta_dist `full_join` and fatty_acid_ppcp_meta_dist `inner_join`
  + Removed the comment `# add the site labels` from the nmds `ggplot` call. Seemed like it was referring to an old version of the code, but I may have misunderstood.
  + Changed `mean_var` definition to include naming inside the function call so avoid the next two lines of renaming
  + I translated the code chunk for all data_scores into mutates
  + Filled in some function arguments
  +
