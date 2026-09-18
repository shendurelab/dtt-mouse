"""Shared core library for the bulk TAPE (DNA Typewriter) analysis pipeline.

Modules:
  io          -- fastq iteration; load/save of references and pickles
  barcode     -- TAPE-BC extraction, Hamming distance, whitelist correction
  parser      -- monomer-junction grammar + robust 5'-forward / 3'-backward array parser
  reference   -- PASS A: derive the TAPE-BC whitelist and insertion vocabulary
  denoise     -- error-collapse + chimera filter + read floor + cross-barcode filter
  depth       -- read-depth model (lambda) -> genomic equivalents / cells
  copynumber  -- per-barcode genomic copy number -> total number of integrations
  editchain   -- the edit-chain trie: the lineage model (prefix tree over edit symbols)
  rarefaction -- Hurlbert interpolation + Chao1 extrapolation of lineage diversity
  tree        -- trie layout + per-integration lineage-tree rendering

The canonical de-noising logic lives ONCE, in denoise.py, and the lineage model
lives ONCE, in editchain.py -- every stage and figure imports them from there.
Both modules, and parser.py, are shared verbatim with the single-cell pipeline:
the read body downstream of the integration barcode is identical in both assays,
so only the barcode read-off differs.
"""
