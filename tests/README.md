# **plant-food-research-open/genepal** Tests

## Minimal Testing

If [Nextflow](https://www.nextflow.io/docs/latest/install.html#install-nextflow) and [Docker](https://docs.docker.com/install) are installed on the system, the pipeline can be minimally tested with the following command:

```bash
nextflow run plant-food-research-open/genepal -r main -profile docker,test --outdir results
```

Or using [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html):

```bash
nextflow run plant-food-research-open/genepal -r main -profile singularity,test --outdir results
```

## nf-test and Continuous Integration (CI)

The GitHub [CI action](../.github/workflows/ci.yml) included with the pipeline continuously tests the pipeline and its components using [nf-test](https://www.nf-test.com). Many components included with the pipeline such as [star/align](../modules/nf-core/star/align) include their own [tests](../modules/nf-core/star/align/tests/main.nf.test) with test data from nf-core. Currently, the full pipeline-level test is run with empty data files in [stub](https://www.nextflow.io/docs/stable/process.html#stub) mode. In this mode data is not processed, rather, the focus is on testing the integrity of data flow through various workflows of the pipeline.

## Testing with a Large Dataset at Plant&Food

Before each release, the functionality of the entire pipeline is tested with a large dataset on the on-prem SLURM-based HPC at The New Zealand Institute of Plant and Food Research.
