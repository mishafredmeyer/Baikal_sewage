### Code cleaning notes and process MRB beginning 2020-02-11

#### 01_data_cleaning.R
**Doing**:

**To do**:
+ Double check old version of piped code output against new version
  (`all_equal`)


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
+ Changed 'PPCP.SUM' to 'ppcp_sum'

#### 02_data_cleaning.R
**Doing**:

**To do**:
+ Should we use `stringsAsFactors = F` for all `read.csv` here and in other scripts
  + MFM: It's up to you. I'll defer to your expertise on that front.

**Done**:
+ Remove commented out code
+ New lines at 80 characters where reasonable
+ Add funtion arguments
+ Comment
+ Line breaks?
+ Can things be purrr'ed?
  + Maybe, but for now I don't think it's worth the hoops you'd have to jump through to do it.
  + MFM: Agreed. I had the same thought, but I am not sure that it is worth the rearrangements.
+ Lint

#### 03_community_composition_analysis.R
**Doing**:

**To do**:
+ Remove commented out code
+ New lines at 80 characters where reasonable
+ Comment
  - Low/high/moderate for sites - what does this mean?
  - MFM: I have added a clarifying comment.
+ Line breaks?
+ Are PI_Groups and POP_GROUPS the same thing?
  - MFM: I have standardized this to read "IDW_pop_group"
+ Standardize the way that "LI-1" is removed from each dataset (i.e., make it happen once per analysis)
  -MFM: This no longer happens. LI-1 remains throughout all analyses.
+ The name data_scores (and species_scores?) is reused - tailor it to each analysis?
  - MFM: I am indifferent. I will defer to your expertise.
+ Could probably improve naming of the several full_joins that happen. Similar names, not descripive enough to differentiate.
  - MFM: If you have an idea of how to do this, I can get behind it.
+ Double check old version of code against new...making some changes that shouldn't affect anything but nervous
  - MFM: Things appear to be working fine.
+ Should section 4 have more subheads?
  - MFM: I have included subheadings.

**Done**:
+ Add funtion arguments
+ Lint
