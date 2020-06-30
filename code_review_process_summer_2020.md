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

#### 02_data_cleaning.R
**Doing**:

**To do**:


**Done**:
 + Edited some formatting (spacing)
 + Changed to include `"stringsAsFactors" = FALSE` in the `read.csv()` statements. Wasn't a huge issue, but I get nervous when joins start leaving messages that factor levels aren't matching and are being automatically reformatted into strings.

#### 03_community_composition_analysis.R
**Doing**:

**To do**:

**Done**:
