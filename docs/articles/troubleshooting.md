# Troubleshooting

**When installing `spatialGE`, I get
`Error in utils::download.file(url, path, method = method, quiet = quiet`.
How can I solve this error?**

This error might appear when installing `spatialGE`:

    Downloading GitHub repo FridleyLab/spatialGE@oospina_dev
    Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
      download from 'https://api.github.com/repos/FridleyLab/spatialGE/tarball/oospina_dev' failed

A workaround is to set `options(timeout)` to a large number:

    options(timeout=9999999) # To avoid R closing connection with GitHub
    devtools::install_github("fridleylab/spatialGE")

**The `gene_interpolation` function outputs an error:
`Error in solve.default(qr.R(qr.VT))`. What can be done about it?**

When running the `gene_interpolation` function with `REML=TRUE`, you may
get this error:

    Error in solve.default(qr.R(qr.VT)) : 
      Lapack routine dgesv: system is exactly singular: U[1,1] = 0
    Error in call to optim
    spatialProcess: Problems with optim in mKrigMLEJoint 
    returned object includes the likelihood evaluations up to the error

This error occurs during REML parameter estimation, and it means that
the gene expression surface estimated might not be accurate.
Nonetheless, if you run the `STplot_interpolation` on the STlist
resulting from `gene_interpolation`, you will notice that gene
expression surfaces have actually been generated for the other genes and
samples despite the error being presented.

If you still want to obtain gene expression surfaces for the other genes
for which `gene_interpolation` could not produce a surface, you can try
`REML=F`, which is a slower procedure.

**When using readRDS to load an STList, I get ‘Error: vector memory
exhausted’. What can I do?**

This error results from R’s default memoery limits. Provided your
computer has enough memory to hold an STList in memory, this error
message can be overridden by typing in a Terminal (Mac OS):

    echo 'R_MAX_VSIZE=100Gb' >> ~/.Renviron

Then, close and re-open R and try to load the STList again.
