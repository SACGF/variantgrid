#!/bin/bash

service gunicorn stop

service celeryd_analysis_workers stop
service celeryd_annotation_workers stop
service celeryd_db_workers stop
service celeryd_web_workers stop
service celeryd_variant_id_single_worker stop
service celeryd_scheduling_single_worker stop
service celeryd_seqauto_single_worker stop
service celeryd_beat stop

pkill gunicorn
