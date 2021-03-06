## 0.0.8 (6 June 2016)

- Clean up BED file prioritization output to be only a list of identified genes.
- Correctly prioritize genes lacking ID fields (nil versus "." issues).

## 0.0.7 (8 April 2016)

- Handle longer BED inputs like BED 6 by truncating to BED 4. Thanks to
  Sven-Eric Schelhorn.
- Do not be as restrictive in filtering larger events on first pass speed up
  since we subset those later and will miss larger ones in ends.

## 0.0.6 (20 March 2016)

- Avoid memory problems with large inputs by avoiding loading
  intersection file into memory.
- Improve speed on large VCF inputs by pre-selecting with bcftools and
  structural variant length.
- Raise Java exceptions like OutOfMemoryError to avoid silent failures.

## 0.0.5 (22 November 2015)

- Carry along both ends of breakends when one falls in a prioritized region.
  Enables resolution of large or inter-chromosomal events originating in a
  region of interest.

## 0.0.4 (20 November 2015)

- Skip reporting of large structural variants with no evidence in breakends.

## 0.0.3 (19 November 2015)

- For large structural variant events (> 50kb) only report genes intersecting
  with end points. Enables more useful annotation for long range events like
  gene fusions.
- Enable build 38 coordinates for preparing known gene BED file from CIViC.

## 0.0.2 (19 October 2015)

- Enable known prioritization with a pre-prepared BED file.
