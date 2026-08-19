# Miniature synthetic reference

This reference is an engineering fixture, not a biological genome. It is generated with a
fixed seed and contains three annotated genes. Its purpose is to make alignment/counting
expectations exact and to permit a very small STAR index for smoke tests.

For STAR genome generation, use `--genomeSAindexNbases 4` and the test configuration's
`sjdb_overhang`. The production GDC reference/index remains unchanged.
