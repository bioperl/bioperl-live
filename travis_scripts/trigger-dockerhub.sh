#!/bin/bash

set -e

# the Docker Hub API endpoint
dockerapi="https://registry.hub.docker.com/u/bioperl/bioperl/trigger"

## Travis runs a build for several versions of Perl, but we want to
## trigger only for the same version of Perl as is running in the
## Docker container
docker_perl="5.18"

if [[ ${TRAVIS_PERL_VERSION:0:4} != "$docker_perl" ]] ; then
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

if [[ -n "$COVERAGE" && "$COVERAGE" != "0" ]] ; then
    echo "Not triggering Docker Hub for code coverage build."
    exit 0
fi

## Should check for tag names that indicate release candidates rather
## than release names, and then skip those.
## However, this is already taken care of by the regular expression
## for whitelisting branches.

if [[ -n "$TRAVIS_TAG" && "$TRAVIS_TAG" != "false" ]] ; then
    # if we are building a whitelisted tag, we trigger the stable image
    echo "Triggering rebuild of Docker image bioperl/bioperl:stable"
    curl -H "Content-Type: application/json" \
         --data '{"docker_tag": "stable"}' \
         -X POST $dockerapi/$DOCKERHUB_TOKEN/
elif [[ "$TRAVIS_BRANCH" != "master" ]] ; then
    # if someone were to create a branch that matches the whitelisting
    # pattern, we skip that here, and only trigger on master
    echo "Not triggering Docker Hub for branches other than master"
else
    # not a pull request, not a tag, and the branch is master
    echo "Triggering rebuild of Docker image bioperl/bioperl:latest"
    curl -H "Content-Type: application/json" \
         --data '{"docker_tag": "latest"}' \
         -X POST $dockerapi/$DOCKERHUB_TOKEN/
fi
