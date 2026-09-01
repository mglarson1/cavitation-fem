# Contributing

Please keep changes reproducible and scoped to the numerical methods in the
accompanying manuscript.

Before proposing a change:

1. run `reproduce('quick')` in MATLAB
2. run the affected table or figure driver directly
3. explain any changed reference value and its numerical cause
4. preserve the documented sampling and cavity-classification protocols

Before a release, also run `reproduce('certify')`. This is intentionally
expensive because it checks the full table matrix, prose-only quantitative
claims, and the safeguarded Newton comparison.

Do not update a reference CSV merely to silence a regression failure. A
changed value should be supported by a corrected method, solver, or extraction
rule and reconciled with the manuscript.
