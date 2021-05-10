# Creates object class.
setClass("STList", slots=list(counts="list",
                              coords="list",
                              clinical="tbl",
                              voom_counts="list",
                              gene_stdev="list",
                              gene_het="list",
                              gene_krige="list",
                              cell_deconv="list",
                              #cell_stdev="list",
                              cell_het="list",
                              cell_krige="list",
                              #purity_krige="list",
                              #prediction_grid="list",
                              prediction_border="list"
),
)
