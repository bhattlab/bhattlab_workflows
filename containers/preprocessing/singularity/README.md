# bhattlab preprocessing singularity container

A singularity container is required to run on SCG. This directory contains
the ```.def``` file and instructions needed to build the singularity container
for common preprocessing tasks. See ../docker/README.md for details of the
contents of the container.

## Building the container

In order to build the container singularity's 'remote' build facility must
be used since local building is not supported on SCG (it requires root privilige,
whereas remote building does not). As its name suggests remote building requires
the use of a remote service, which for singularity is Sylabs.io. Sylabs
are the developers and maintainers of singularity and offer a build service
and a container registry which includes a generous free tier. In order to
use the service you need to create an account on Sylabs.io and to generate
an access token. You can use your github account to login to Sylabs.io and
you can create an access token [here](https://cloud.sylabs.io/dashboard#security).
Once you have an access token, create a configuration file in ```~/.singularity/remote.yaml```

```yaml
Active: SylabsCloud
Remotes:
  SylabsCloud:
    URI: cloud.sylabs.io
    Token: <your access token>
    System: true
    Exclusive: false
```

The definition file ```preprocessing-docker-micromamba.def``` is akin to
the Dockerfile used to build the docker container, but for singularity. Note
that the ```From:``` statement must refer to the same container built and pushed
from the docker directory (../docker/docker-image-name).

Given a remote.yaml configuration file and the definition file, the container
can be built with the following command:

```singularity build --remote bhattlab-micromamba-preproc.img preprocessing-docker-micromamba.def```

