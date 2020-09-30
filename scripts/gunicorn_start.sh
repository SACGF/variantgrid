#!/bin/bash

PORT=${PORT:-8000}
#IP_ADDRESS=$(ifconfig eth0 | grep 'inet addr:' | cut -d: -f2 | awk '{ print $1}')
IP_ADDRESS=${IP_ADDRESS:-localhost}:${PORT}

cd "$(dirname $0)/.."
gunicorn --bind ${IP_ADDRESS} -t 300 -w 4 --log-level=DEBUG variantgrid.wsgi:application

