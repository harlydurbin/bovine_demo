# Demographic and selection history of cattle and related species

## Meta-data processing

1. SRA metadata scraped from downloaded XML files then initially tidied in `notebooks/xml_scraping.Rmd`
2. Population labels futher categorized by region, continent, and species in `notebooks/categorizing.Rmd`
3. Coverage of samples in tidied dataset assessed `notebooks/coverage.Rmd`
4. Duplicate samples identified using `source_functions/filter_eval.snakefile` then removed from metadata in `source_functions/duplicate_samples.R`

**Resulting sample metadata file: `data/derived_data/sample_metadata.csv`**