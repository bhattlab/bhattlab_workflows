# Docker container bhattlab preprocessing

This directory contains the Dockerfile and associated files for the bhattlab
preprocessing container. It uses micromamba rather than conda to create the
preprocessing environment to minimized build time and image size. The
resulting image is about half the size of the equivalent conda image
and builds in a minute or so vs 10s of minutes.

The container is pushed to the bhattlab github container registry (called
packages on github). This docker container can then be used to build
a singularity container on SCG (see ../singularity) for details.

Building the container is straightforward on any machine with docker installed
via the ```build.sh``` script. In order to push the container to the
github container registry a github token is required
([from here](https://github.com/settings/tokens)). Note that a 
"personal access token '(cliassic)'" is required rather than the newer granular tokens.
This token can be used to login to the github registry as follows:

```bash
echo "<token>" | docker login ghcr.io -u "<your-github-username>" --password-stdin
```

The script ```push.sh``` can then be used to push the container to github.

The package can be viewed on github [here](https://github.com/orgs/bhattlab/packages?repo_name=bhattlab_workflows). The visibility of a new 'package' must be
changed manually to public since the default is private.

