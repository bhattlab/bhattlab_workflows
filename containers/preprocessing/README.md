# bhattlab preprocessing container

This directory contains the scripts needed to build the bhattlab preprocessing container.
The container is built in two steps:
1. A docker container is built that contains the desired tools. Docker is preferred
for this step since it is both faster and easier to run on a laptop than
singularity. This container is pushed to GitHub's container registry
(GitHub Packages) to make it accessible from SCG. This step must be run
on a machine with docker installed, typically your laptop.
2. The docker container is used to build a singularity container for use
on SCG. This step is typically run on an SCG login node.

To run these tools, accounts on GitHib and Sylabs.io (developers/maintainers
of singularity) are required; as are access tokens from each of these.
Accounts and tokens are free, or have generous free tiers.

Each step is described in more detail in the sudirectories for docker and
singularity, please see the README.md files in those directories for more details.

See docker/Dockerfile for the contents of the container.