# bcbio.prioritize

Prioritize small variants, structural variants and coverage based on biological
inputs. The goal is to use pre-existing knowledge of relevant genes, domains and
pathways involved with a disease to extract the most interesting signal from a
set of high quality small or structural variant calls. Given information on
coverage, it will be able to identify poorly covered regions in potential genes
of interest.

This is exploratory work in progress.

## Usage

Create a database of priority regions based on gene and domain regions with
biological evidence:

    bcbio-prioritize createdb known.db transcript.bed known_input.vcf known_input.bed

Prioritize a set of structural variant calls in BED format given a database of
regions:

    bcbio-prioritize known known.db calls.bed > calls-known.bed

Identify regions missing sufficient sequencing coverage:

    bcbio-prioritize missing known.db coverage.db > known-missing.bed

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
