#!/bin/bash

cd "$(dirname $0)/.."
IP_ADDRESS=$(ifconfig eth0 | grep 'inet addr' | cut -d: -f2 | awk '{print $1}')
gunicorn --bind ${IP_ADDRESS} -t 300 -w 4 --log-level=DEBUG variantgrid.wsgi:application

