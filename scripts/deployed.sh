#!/bin/bash
ACCESS_TOKEN=2c82541df93e485e9a8ed25160abdf38
ENVIRONMENT=${HOSTNAME//-/}
ENVIRONMENT=${ENVIRONMENT,,}
LOCAL_USERNAME=`whoami`
REVISION=`git rev-parse --verify HEAD`
echo curl https://api.rollbar.com/api/1/deploy/ \
  -F access_token=$ACCESS_TOKEN \
  -F environment=$ENVIRONMENT \
  -F revision=$REVISION \
  -F local_username=$LOCAL_USERNAME

curl https://api.rollbar.com/api/1/deploy/ \
  -F access_token=$ACCESS_TOKEN \
  -F environment=$ENVIRONMENT \
  -F revision=$REVISION \
  -F local_username=$LOCAL_USERNAME

python3 manage.py deployed