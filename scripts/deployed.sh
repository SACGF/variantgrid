#!/bin/bash
ACCESS_TOKEN=908d7bc42049497494d408f7fb9ee6a5
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

python3.8 manage.py deployed