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

## Known gaps as of this writing

- `parallel::detectCores()` returns `NA` in network-sandboxed shells
  (`sysctl` is blocked), so `count_cores()` returns `NA` and any call
  relying on the default `cores = NULL` hard-errors. Pass `cores =`
  explicitly in that situation.
- `src/` will not compile without `RcppEigen` and `RcppProgress` present;
  `install_deps.R` installs both.
