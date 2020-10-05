# Add files 50MB or larger to .gitignore (https://stackoverflow.com/questions/4035779/gitignore-by-file-size)
find . -size +50M | sed 's|^\./||g' | cat >> .gitignore
