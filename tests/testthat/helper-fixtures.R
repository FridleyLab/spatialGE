## Sourced automatically by testthat before running this package's tests.
## Materializes the synthetic fixtures (inst/testdata/make_fixtures.R) into a
## session tempdir once; individual test files just reference `fixtures$...`.

source(system.file("testdata", "make_fixtures.R", package = "spatialGE"), local = TRUE)

.fixture_root <- file.path(tempdir(), "spatialge-fixtures")
dir.create(.fixture_root, showWarnings = FALSE, recursive = TRUE)
fixtures <- fixture_materialize_all(.fixture_root)
