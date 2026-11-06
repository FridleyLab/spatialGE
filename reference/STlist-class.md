# Definition of an STlist object class.

Definition of an STlist object class.

## Slots

- `counts`:

  per spot RNA counts

- `spatial_meta`:

  per spot x,y coordinates

- `gene_meta`:

  per gene statistics (e.g., average expression, variance, Moran's I)

- `sample_meta`:

  dataframe with metadata per sample

- `tr_counts`:

  transfromed per spot counts

- `gene_krige`:

  results from kriging on gene expression

- `misc`:

  Parameters and images from ST data
