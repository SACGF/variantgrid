#!/bin/bash

RABBIT_MQ_URL="amqp://guest:guest@localhost"
REDIS_URL="redis://localhost:6379/0"
BROKER_URL=${RABBIT_MQ_URL}

cd "$(dirname $0)/.." || exit;

# One CPU-bound VEP per annotation slot: keep this x ANNOTATION_VEP_FORK at or under the box's cores (#2477)
ANNOTATION_WORKERS_CONCURRENCY=${ANNOTATION_WORKERS_CONCURRENCY:-4}

celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=4 -n analysis_workers --queues analysis_workers &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=${ANNOTATION_WORKERS_CONCURRENCY} -n annotation_workers --queues annotation_workers &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=4 -n db_workers --queues db_workers &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=4 -n web_workers --queues web_workers &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=1 -n variant_id_single_worker --queues variant_id_single_worker &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=1 -n scheduling_single_worker --queues scheduling_single_worker &
celery -b ${BROKER_URL} --app variantgrid worker -l DEBUG --max-tasks-per-child=100 --concurrency=2 -n heavy_workers --queues heavy_workers &
wait
