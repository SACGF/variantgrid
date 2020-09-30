#!/bin/bash

cd "$(dirname $0)/.." || exit;

# more experimentation with event queues and concurrency is needed
# everything is currently set to be single threaded
RABBIT_MQ_URL="amqp://guest:guest@localhost"
REDIS_URL="redis://localhost:6379/0"
BROKER_URL=${RABBIT_MQ_URL}

celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P gevent -n analysis_workers --queues analysis_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P gevent -n annotation_workers --queues annotation_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P gevent -n db_workers --queues db_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P gevent -n web_workers --queues web_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P gevent -n variant_id_single_worker --queues variant_id_single_worker &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P eventlet -n scheduling_single_worker --queues scheduling_single_worker &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -P eventlet -n seqauto_single_worker --queues seqauto_single_worker &
wait
