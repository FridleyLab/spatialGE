# Development environment

spatialGE development happens inside the `spatialge-dev` conda environment,
not system R, because the geospatial/Bioconductor stack (`sf`, `gdal`,
`hdf5r`, `EBImage`, ...) is easier to satisfy from conda-forge than from
CRAN source builds.

## Activate

```sh
conda activate spatialge-dev
export R_MAKEVARS_USER="$(pwd)/dev/environment/Makevars"
```

`spatialGE-dev.yml` is the conda spec used to (re)create the environment:

```sh
conda env create -n spatialge-dev -f dev/environment/spatialGE-dev.yml
```

After the environment exists, install the R package dependencies (needs
network):

```sh
Rscript dev/environment/install_deps.R
```

This records the resulting toolchain in `dev/environment/RECORDED.md`
(`sessionInfo()`, `extSoftVersion()`, BLAS/LAPACK).

## renv

`renv/`, `renv.lock`, and `.Rprofile` are all gitignored. `.Rprofile`
would otherwise source `renv/activate.R` on every R startup; on a fresh
clone with no `renv/library/`, that autoloader reaches for the network and
hangs under a sandboxed/offline shell. Reproducibility here comes from the
conda yaml plus `install_deps.R`, not from an renv lockfile.

If you deliberately re-enable renv for local work and hit the same hang,
start R with the autoloader off:

```sh
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE R
```

## Corporate network / TLS interception

On networks behind an SSL-inspecting proxy (e.g. Netskope), `curl`/R's
`download.file()`/`conda` will fail with "SSL certificate problem:
self-signed certificate in certificate chain" (curl exit code 60) because
the proxy's root CA isn't in the default trust store. Fix by building a
combined bundle of the system roots plus the corporate root(s) pulled from
the macOS System keychain, then pointing every tool at it:

```sh
security find-certificate -c "<corporate root CN>" -p \
  /Library/Keychains/System.keychain >> corp_roots.pem
cat /etc/ssl/cert.pem corp_roots.pem > combined_cacert.pem

export CURL_CA_BUNDLE=/path/to/combined_cacert.pem
export SSL_CERT_FILE=/path/to/combined_cacert.pem
export REQUESTS_CA_BUNDLE=/path/to/combined_cacert.pem
```

`CURL_CA_BUNDLE`/`SSL_CERT_FILE` cover `curl` and R's own downloaders;
`REQUESTS_CA_BUNDLE` covers conda. Export all three before running
`conda env create`/`conda install` or `install_deps.R`.

### Concrete values observed on the Moffitt network

The intercepting root is `O=Moffitt Cancer Center, CN=ca.moffitt.goskope.com`,
so the first command above is:

```sh
security find-certificate -a -c ca.moffitt.goskope.com -p \
  /Library/Keychains/System.keychain > corp_roots.pem
```

Note that macOS trusts this CA via the system keychain, so `/usr/bin/curl`
succeeds while conda's `curl` and conda R's libcurl both fail -- "curl works, so
the network is fine" is a false conclusion here. `CURL_CA_BUNDLE=/etc/ssl/cert.pem`
alone is also insufficient; the corporate root has to be appended. If the System
keychain lookup comes back empty, try the login keychain
(`~/Library/Keychains/login.keychain-db`) as well.

Verify from conda R itself before trusting any network-dependent run. Expect
`HTTP/2 302`, not `libcurl error code 60`:

```sh
Rscript -e 'cat(curlGetHeaders("https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=")[1])'
```

### Why this is worth checking explicitly

Left unfixed, the failure is close to silent and is actively misattributed
rather than loud:

- `tests/testthat/test-tests.R` calls `skip_if_offline()`, so the suite passes
  having quietly run one test fewer.
- `vignettes/basic_functions_vignette.Rmd` guards ~30 chunks on `eval=net_check`,
  so a failed download renders a green but essentially empty document.
- Every `@examples` block reports "Could not run example. Are you connected to
  the internet?" regardless of the real cause -- that same message also appears
  for genuine code defects, which is exactly how the `zip_tmp` bug fixed in
  2.0.0 stayed hidden across two releases.

So a green run with a broken CA bundle proves much less than it appears to.

## Known gaps as of this writing

- `parallel::detectCores()` returns `NA` in network-sandboxed shells
  (`sysctl` is blocked), so `count_cores()` returns `NA` and any call
  relying on the default `cores = NULL` hard-errors. Pass `cores =`
  explicitly in that situation.
- `src/` will not compile without `RcppEigen` and `RcppProgress` present;
  `install_deps.R` installs both.
