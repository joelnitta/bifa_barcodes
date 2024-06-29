# bifa_barcodes

Code repository for ["Improving knowledge of Asian pteridophytes through DNA sampling of specimens in regional collections" project](https://www.gbif.org/project/BIFA6_010/improving-knowledge-of-asian-pteridophytes-through-dna-sampling-of-specimens-in-regional-collections)

There are two options (Docker images) provided to run the code: either [run everything from scratch](#docker-image-to-run-the-analysis-from-scratch), or launch an RStudio session with the [already completed workflow](#docker-image-with-completed-analysis).

## Docker image to run the analysis from scratch

First, clone this GitHub repo then enter it:

```
git clone https://github.com/joelnitta/bifa_barcodes
cd bifa_barcodes
```

To run the analysis pipeline, do:

```
docker run --rm -dt -v ${PWD}:/wd -w /wd joelnitta/bifa_barcodes:latest Rscript -e "targets::tar_make()"
```

### Interacting with the Docker container while it runs

To inspect the output of each workflow step, first launch a Docker container in the background:

```
docker run --rm -dt -v ${PWD}:/wd -w /wd joelnitta/bifa_barcodes:latest
```

Then obtain the Docker container ID `docker ps` and enter it with `docker exec -it <CONTAINER ID>`.

## Docker image with completed analysis

Another Docker image is provided that already contains the complete workflow *after* everything has finished running. That way, you don't have to wait for each step to finish. Also, this one comes with RStudio, so you can explore the output in RStudio in your browser.

### Start up RStudio

First, run this command in the terminal:

```
docker run --rm --name bifa -dt -e PASSWORD=pw -p 8787:8787 joelnitta/bifa_barcodes_demo:latest
```

Next, navigate to <http://localhost:8787/> in your web browser.

In the login screen, enter `rstudio` for user and `pw` for password.

This should take you to an instance of RStudio in your browser containing the complete code base, raw data, and results.

### Browse the code

I recommend first running `source("_targets.R")` to load all the necesssary packages.

Open the `_targets.R` file to see the analysis pipeline. The results (targets) of individual steps can be loaded into memory with `tar_load()`, for example: `tar_load(ppgi)` loads the dataframe with the PPG I classification of pteridophytes.

The barcoding analysis section of the manuscript was produced with the Quarto markdown file in `reports/index.qmd`. You can also open that file and step through the code chunks to see how the tables and graphs were created.

### Close down the container

When you are done, go back to the terminal and run `docker kill bifa` to close down the container.