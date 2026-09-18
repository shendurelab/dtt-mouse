"""Shared library for the bulk TAPE analysis pipeline.

Modules:
  io        -- fastq iteration; load/save of references and pickles
  barcode   -- TAPE-BC extraction, Hamming distance, whitelist correction
  parser    -- monomer-junction classification + robust 3'-recovery array parser
  reference -- PASS A: derive the TAPE-BC whitelist and NNN insertion vocabulary
  denoise   -- error-collapse + UCHIME chimera filter + cross-barcode filter
  depth     -- read-depth model (lambda) -> genomic equivalents / cells
  tree      -- prefix-trie construction over ordered edits (per-integration trees)

The canonical de-noising logic lives ONCE, in denoise.py -- every stage and
figure imports it from there (the original workspace had four divergent copies).
"""
