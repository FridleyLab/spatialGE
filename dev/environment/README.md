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

**There is more than one interceptor, and which one you get depends on where you
are.** An earlier version of this section named only
`CN=ca.moffitt.goskope.com` and told you to extract that one cert. That is not
enough: it is an *intermediate*, not a root, and a different network path
intercepts with an entirely different CA. Two chains observed, both on Moffitt
networks:

| Where | Chain presented | Trust anchor needed |
|---|---|---|
| Netskope cloud proxy | `<host>` → `ca.moffitt.goskope.com` → `caadmin.netskope.com` | `caadmin.netskope.com` |
| On-premises | `<host>` → `Moffitt Issuing Authority` → `Moffitt Root CA` | `Moffitt Root CA` |

Interception is also **per host, not per network**: on the cloud-proxy path
`conda.anaconda.org` was intercepted while `github.com` was not, so R's
`download.file()` from GitHub worked while `conda search` failed. On the
on-premises path `doi.org` and `satijalab.org` are intercepted, which makes
`R CMD check --as-cran` report every DOI citation as a possibly invalid URL
(with `libcurl error code 60` and, tellingly, `(Status without verification:
OK)`).

So do not extract one CN by name. Take every corporate cert the keychain holds,
roots and intermediates, and let the extras be harmless:

```sh
security find-certificate -a -p /Library/Keychains/System.keychain > sys.pem

python3 - <<'EOF' > corp_roots.pem
import re, subprocess
pem = open('sys.pem').read()
for c in re.findall(r'-----BEGIN CERTIFICATE-----.*?-----END CERTIFICATE-----', pem, re.S):
    subj = subprocess.run(['openssl', 'x509', '-noout', '-subject'],
                          input=c, capture_output=True, text=True).stdout.lower()
    if 'moffitt' in subj or 'netskope' in subj:
        print(c)
EOF

cat /etc/ssl/cert.pem corp_roots.pem > combined_cacert.pem
```

If you prefer to do the filtering by hand, note that `security
find-certificate -c "Moffitt Root CA"` returns **nothing** even though that
cert is present -- `-c` does not match its subject. Enumerate with `-a` and
filter on the decoded subject instead. On this machine that yields 11 corporate
certs, taking the bundle from 129 to 139.

Verify the result is a superset before replacing a working bundle, by comparing
SHA-256 fingerprints of the old and new files; losing a public root here breaks
far more than it fixes.

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
