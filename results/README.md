# Results policy

`results/r1/` is reserved for fresh major-revision receipts. Historical submitted-paper numbers live in `validation/reference_values.json` and are never treated as fresh evidence.

A receipt should contain at minimum:

- gate id;
- UTC timestamp;
- host/OS/Python identity;
- Git commit;
- scientific fingerprint;
- input hashes;
- observed values;
- reference values and tolerance;
- PASS/FAIL;
- command and exit code.
