#!/bin/bash

# Runs a command line under perlbrew in home dir
source ${HOME}/perl5/perlbrew/etc/bashrc 
perlbrew switch 5.38.0
exec "$@"
