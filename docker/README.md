This file contains simple instructions for how you can develop chimeraviz, either the current release verion of the development version, using the RStudio Docker images.

If these instructions don't work for you, please don't hesitate to [open an issue](https://github.com/stianlagstad/chimeraviz/issues/new) and I'll assist you.

# Release version
1. Build the Docker image with `docker build . -f Dockerfile.release -t chimeraviz_release` (this may take several minutes)
2. Start a container with `docker run -it -e PASSWORD=42 -p 8787:8787 -v /path/to/chimeraviz/:/chimeraviz chimeraviz_release`. Remember to change `/path/to/chimeraviz/` to where you cloned the chimeraviz repository.
3. Go to localhost:8787, login with rstudio/42, run `setwd("/chimeraviz")`, then `devtools::load_all()`, and start coding

# Development version
1. Build the Docker image with `docker build . -f Dockerfile.devel -t chimeraviz_devel` (this may take several minutes)
2. Start a container with `docker run -it -e PASSWORD=42 -p 8787:8787 -v /path/to/chimeraviz/:/chimeraviz chimeraviz_devel`. Remember to change `/path/to/chimeraviz/` to where you cloned the chimeraviz repository.
3. Go to localhost:8787, login with rstudio/42, run `setwd("/chimeraviz")`, then `devtools::load_all()`, and start coding

