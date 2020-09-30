#!/bin/bash

RABBIT_MQ_URL="amqp://guest:guest@localhost"
REDIS_URL="redis://localhost:6379/0"
BROKER_URL=${RABBIT_MQ_URL}

cd "$(dirname $0)/.." || exit;

celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=4 -n analysis_workers --queues analysis_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=4 -n annotation_workers --queues annotation_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=4 -n db_workers --queues db_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=4 -n web_workers --queues web_workers &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -n variant_id_single_worker --queues variant_id_single_worker &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -n scheduling_single_worker --queues scheduling_single_worker &
celery worker -l DEBUG -b ${BROKER_URL} --app variantgrid --concurrency=1 -n seqauto_single_worker --queues seqauto_single_worker &
wait
