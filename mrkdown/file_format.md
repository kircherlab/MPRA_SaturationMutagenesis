### Format description

Variant files for each element are available for genome releases GRCh37 and GRCh38. Except of the chromosomal coordinates the files of the different releases are identical. Files are TAB (`.tsv`) or COMMA (`.csv`) separated, each row contains one variant, and they include a header. If a possible variant at a position is not shown, it was not observed in the saturation mutagenesis library and therefore not present in the final model for fitting. 

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

### Data usage

We are making our data available prior to publication in line with  [Fort Lauderdale principle](https://www.genome.gov/pages/research/wellcomereport0303.pdf), allowing others to use the data but allowing the data producers to make the first presentations and to publish the first paper with global analyses of the data. In addition, we also reserve the right to publish the first analysis of the differences seen in the TERT knock-down experiments and alternative cell-type experiments. Studies that do not overlap with these intentions may be submitted for publication at any time, but must appropriately cite the data source. After publication of the data, the first publication of the data producers should be cited for any use of this data. 
