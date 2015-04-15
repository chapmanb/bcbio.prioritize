# bcbio.prioritize

Prioritize small variants, structural variants and coverage based on biological
inputs. The goal is to use pre-existing knowledge of relevant genes, domains and
pathways involved with a disease to extract the most interesting signal from a
set of high quality small or structural variant calls. Given information on
coverage, it will be able to identify poorly covered regions in potential genes
of interest.

This is exploratory work in progress.

## Usage

Create file of priority regions based on gene and domain regions with biological
evidence:

    bcbio-prioritize create -o known.bed.gz -b transcript.bed -k input1.vcf -k input2.bed

Prioritize a set of structural variant calls in BED format given a binned set of known
changes:

    bcbio-prioritize known -i calls.bed.gz -k known.bed.gz -o calls-known.bed.gz

Identify regions missing sufficient sequencing coverage:

    bcbio-prioritize missing known.db coverage.db > known-missing.bed

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
