# bifa_barcodes

Code repository for ["Improving knowledge of Asian pteridophytes through DNA sampling of specimens in regional collections" project](https://www.gbif.org/project/BIFA6_010/improving-knowledge-of-asian-pteridophytes-through-dna-sampling-of-specimens-in-regional-collections)

## Docker image

The Docker image `joelnitta/bifa_barcodes:latest` is automatically built and pushed to Docker Hub via a [Github Workflow](.github/workflows/docker.yaml) .

## Running the analysis

To run the analysis pipeline (requires Docker), do:

```
docker run --rm -dt -v ${PWD}:/wd -w /wd joelnitta/bifa_barcodes:latest Rscript -e "targets::tar_make()"
```

## Interacting with the Docker container

To inspect the output of each workflow step, first launch a Docker container in the background:

```
docker run --rm -dt -v ${PWD}:/wd -w /wd joelnitta/bifa_barcodes:latest
```

Then obtain the Docker container ID `docker ps` and enter it with `docker exec -it <CONTAINER ID>`.
