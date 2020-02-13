### Code cleaning notes and process MRB beginning 2020-02-11

#### 01_data_cleaning.R
**Doing**:

**To do**:
+ Double check old version of piped code output against new version (`all_equal`)


**Done**:
+ Add standardized section headers with consistent spacing and caps
+ Add in function arguments
+ Standardize how `unite()` is called
+ Add line breaks between variables in the script
+ There's a `high` vector being referenced in the periphyton section but I can't find any history of it in script 01. Figure out what's supposed to be here. **update**: MM said erroneously included in this script.
+ Commenting
+ Remove dates from filenames
+ Break pipe chains with >= 10 pipes into separate objects
+ Linting

#### 02_data_cleaning.R
**Doing**:

**To do**:
+ Should we use `stringsAsFactors = F` for all `read.csv` here and in other scripts?

**Done**:
+ Remove commented out code
+ New lines at 80 characters where reasonable
+ Add funtion arguments
+ Comment
+ Line breaks?
+ Can things be purrr'ed?
  + Maybe, but for now I don't think it's worth the hoops you'd have to jump through to do it.
+ Lint

#### 03_community_composition_analysis.R
**Doing**:

**To do**:
+ Remove commented out code
+ New lines at 80 characters where reasonable
+ Add funtion arguments
+ Comment
+ Line breaks?
+ Can things be purrr'ed?
  + Maybe, but for now I don't think it's worth the hoops you'd have to jump through to do it.
+ Lint
+
**Done**:
