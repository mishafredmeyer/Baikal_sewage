### Code review notes and process MRB

#### 01_data_cleaning.R
**Doing**:
+ New lines at 80 characters where possible
+ Standardize how `unite()` is called
+ Add spacing between variables in the script

**To do**:
+ Remove dates from filenames
+ Break pipe chains with >= 10 pipes into separate objects
+ Check old version of piped code output against new version (all_equal)

**Done**:
+ Add standardized section headers with consistent spacing and caps
+ Add in function arguments
