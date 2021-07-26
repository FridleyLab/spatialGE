##
# Definition of an STList object class.
#
#
setClass("STList", slots=list(counts="list",
                              coords="list",
                              clinical="tbl",
                              voom_counts="list",
                              gene_stdev="list",
                              gene_het="list",
                              gene_krige="list",
                              cell_deconv="list",
                              cell_het="list",
                              cell_krige="list",
                              gene_krige_data="list",
                              deconv_krige_data="list",
                              st_clusters="list"
),
)
