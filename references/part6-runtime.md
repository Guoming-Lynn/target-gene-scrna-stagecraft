# Part 6 runtime

Use an isolated environment matching the frozen official source revision and model.
The older dependency list in part6-requirements.txt is a historical candidate,
not a verified installation recipe or an official compatibility guarantee.
Record resolved versions after running imports, model load and official parity.

Reference source revision: 04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5.
Reference model: Geneformer-V2-104M, maximum sequence length 4096.
Freeze hashes of model weights, configuration and tokenizer dictionaries.
Set STAGECRAFT_GENEFORMER_SOURCE and STAGECRAFT_GENEFORMER_MODEL for the probe;
install the source in the execution environment or set PYTHONPATH for subsequent commands.

KO requires target presence in the final official token sequence. Record truncated
cells separately. OE comparator overflow must match the pinned official helper.
Preserve per-cell parity ledgers. Results are exploratory embedding shifts only;
they do not predict expression or establish a causal effect.
