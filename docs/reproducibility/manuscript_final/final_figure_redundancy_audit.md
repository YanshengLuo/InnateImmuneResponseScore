# Final figure redundancy audit

Overall status: **PASS**. The four retained views share validated Delta IMRSz sources where appropriate but differ in scope, statistical unit, axes, encoding, and scientific purpose.

| Panel | Data scope | Statistical unit | X variable | Y variable | Graphical encoding | Scientific question | Reason retained |
|---|---|---|---|---|---|---|---|
| Figure 1B | All 68 valid scored contrasts, compacted to the full dataset/context landscape | Dataset/context mean | Mean Delta IMRSz on pseudo-log scale | Compact dataset/tissue/context label | Sized and group-colored point landscape | What does the entire IMRS response landscape look like? | Global orientation across every manuscript group. |
| Figure 3B | Fourteen nested splits from GSE119119, GSE139529, and GSE279743 only | Nested split within independent dataset, plus dataset mean | Split-level Delta IMRSz | Three independent primary datasets | Split points, observed range line, and larger mean diamond | Do the three independent primary datasets each show the same positive transfer pattern? | Makes independence and within-dataset consistency explicit. |
| Figure 4A | Thirteen context-shifted contrasts only | Context-shifted contrast | Observed Delta IMRSz | Biological context category | Jittered manuscript-group-colored points with selected provenance labels | In what biological contexts does the acute interpretation weaken? | Shows biological boundaries rather than another dataset forest. |
| Supplementary Figure S1A | All valid dataset/tissue/time summaries across four manuscript groups | Dataset/context mean at explicit time | Mean Delta IMRSz on pseudo-log scale | Explicit dataset/tissue/time label within manuscript-group facets | Detailed faceted point audit | What is the detailed dataset/context provenance across all groups? | Provides provenance detail beyond the compact main-text overview. |

No pair of main-text panels uses the same rows, axes, and graphical encoding with only a subset change.
