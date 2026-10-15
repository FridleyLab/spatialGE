##
# load_images(): image-to-sample wiring.
#
# Regression coverage for the previously-untested paths. `load_images()` had no
# committed test references at all -- the only exercise was in the gitignored,
# CI-not-run golden tier, which passed a *single directory* whose files were
# named exactly `<sample>.png`. That combination hid two defects:
#
#   1. `if(length(images) == 1 & dir.exists(images))` used the vectorized `&`,
#      so passing more than one path -- which is what the function's own error
#      message tells users to do -- failed with "the condition has length > 1".
#   2. Image decoding was chosen from the file extension only, so a mislabelled
#      image produced an opaque decoder error naming no file.
#
# These tests use the synthetic fixtures, so they are offline and deterministic.

.img_test_stlist <- function() {
  f <- fixtures$list_dfs
  STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, cores = 1)
}

.img_test_samples <- function() fixtures$base$samples

## Writes one small PNG per sample into `dir`, named `<sample>.png`.
.img_test_write <- function(dir, samples = .img_test_samples(), ext = "png") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in samples) fixture_write_image(file.path(dir, paste0(s, ".", ext)))
  dir
}

testthat::test_that("load_images() accepts a single directory", {
  st <- .img_test_stlist()
  d <- .img_test_write(file.path(tempdir(), "li-dir"))

  out <- load_images(st, images = d)

  testthat::expect_s4_class(out, "STlist")
  testthat::expect_setequal(names(out@misc$sp_images), .img_test_samples())
})

testthat::test_that("load_images() accepts a vector of file paths", {
  ## Regression: this used to fail with "the condition has length > 1".
  st <- .img_test_stlist()
  d <- .img_test_write(file.path(tempdir(), "li-vec"))
  paths <- file.path(d, paste0(.img_test_samples(), ".png"))

  testthat::expect_length(paths, length(.img_test_samples()))
  out <- load_images(st, images = paths)

  testthat::expect_setequal(names(out@misc$sp_images), .img_test_samples())
})

testthat::test_that("load_images() accepts several directories in one call", {
  ## The SMI vignette needs images from two per-slide folders; previously that
  ## required one call per folder.
  st <- .img_test_stlist()
  samples <- .img_test_samples()
  d1 <- .img_test_write(file.path(tempdir(), "li-multi-1"), samples[1])
  d2 <- .img_test_write(file.path(tempdir(), "li-multi-2"), samples[-1])

  out <- load_images(st, images = c(d1, d2))

  testthat::expect_setequal(names(out@misc$sp_images), samples)
})

testthat::test_that("load_images() warns, naming the sample, when no image matches", {
  ## e.g. Visium exports that name images after a GEO accession rather than the
  ## sample folder: "GSM6433585_092A_tissue_hires_image.png" vs "sample_092a".
  st <- .img_test_stlist()
  d <- file.path(tempdir(), "li-unmatched")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  fixture_write_image(file.path(d, "GSM0000000_notasample_tissue_hires_image.png"))

  testthat::expect_warning(
    out <- load_images(st, images = d),
    "No image found for sample"
  )
  testthat::expect_length(out@misc$sp_images, 0)
})

testthat::test_that("load_images() honours an explicit sample -> file mapping", {
  ## The escape hatch for filenames that do not resemble the sample names.
  st <- .img_test_stlist()
  samples <- .img_test_samples()
  d <- file.path(tempdir(), "li-named")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  files <- file.path(d, paste0("GSM000", seq_along(samples), "_image.png"))
  for (fp in files) fixture_write_image(fp)

  mapping <- stats::setNames(files, samples)
  out <- load_images(st, images = mapping)

  testthat::expect_setequal(names(out@misc$sp_images), samples)
})

testthat::test_that("load_images() loads an image whose extension is wrong", {
  ## Decoding falls back to the other codec instead of trusting the extension.
  st <- .img_test_stlist()
  d <- .img_test_write(file.path(tempdir(), "li-mislabelled"), ext = "jpg")

  out <- testthat::expect_silent(load_images(st, images = d))

  testthat::expect_setequal(names(out@misc$sp_images), .img_test_samples())
})

testthat::test_that("load_images() errors clearly when given no usable image", {
  st <- .img_test_stlist()

  testthat::expect_error(
    load_images(st, images = file.path(tempdir(), "li-does-not-exist.png")),
    "existing .png, .jpg or .jpeg"
  )
})

testthat::test_that("load_images() still requires an `images` argument", {
  st <- .img_test_stlist()

  testthat::expect_error(load_images(st), "provide a vector with images file paths")
})
