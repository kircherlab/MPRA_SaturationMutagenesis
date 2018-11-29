### Format description

Variant files for each element are available for genome releases GRCh37 and GRCh38. Except of the chromosomal coordinates the files of the different releases are identical. Files are TAB (`.tsv`) or COMMA (`.csv`) separated, each row contains one variant, and they include a header. If a possible variant at a position is not shown it was not present in inserted through error-prone PCR in our target sequence and therefore not present in the final model for fitting. 

#### Columns

**Chromosome** - Chromosome of the variant.  
**Position** - Chromosomal position (GRCh38 or GRCH38) of the variant.  
**Ref** - Reference allele of the variant (`A`, `T`, `G`, or `C`).  
**Alt** - Alternative allele of the variant (`A`, `T`, `G`, or `C`). One base-pair deletions are represented as `-`.  
**Tags** - Number of unique tags associated with the variant.  
**DNA** - Count of DNA sequences that contain the variant (used for fitting the linear model).   
**RNA** - Count of RNA sequences that contain the variant (used for fitting the linear model).  
**Value** - Log2 variant expression effect derived from the fit of the linear model (coefficient).  
**P-Value** - P-value of the coefficient.