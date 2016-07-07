#!/bin/bash

set -e

# the Docker Hub API endpoint
dockerapi="https://registry.hub.docker.com/u/bioperl/bioperl/trigger"

## Travis runs a build for several versions of Perl, but we want to
## trigger only for the same version of Perl as is running in the
## Docker container
docker_perl="5.018"

if [[ ${TRAVIS_PERL_VERSION:0:5} != "$docker_perl" ]] ; then
    echo "Triggering Docker Hub only for Perl $docker_perl, not $TRAVIS_PERL_VERSION"
    exit 0
fi

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]] ; then
    echo "Not triggering Docker Hub for pull requests."
    exit 0
fi

if [[ -z "$DOCKERHUB_TOKEN" ]] ; then
    echo "No API token for Docker Hub, add to repository settings."
    exit 1
fi

## Should check for tag names that indicate release candidates rather
## than release names, and then skip those.

if [[ -n "$TRAVIS_TAG" && "$TRAVIS_TAG" != "false" ]] ; then
    echo "Triggering rebuild of Docker image bioperl/bioperl:stable"
    curl -H "Content-Type: application/json" \
         --data '{"docker_tag": "stable"}' \
         -X POST $dockerapi/$DOCKERHUB_TOKEN/
else
    echo "Triggering rebuild of Docker image bioperl/bioperl:latest"
    curl -H "Content-Type: application/json" \
         --data '{"docker_tag": "latest"}' \
         -X POST $dockerapi/$DOCKERHUB_TOKEN/
fi
