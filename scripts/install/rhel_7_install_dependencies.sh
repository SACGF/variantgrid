#!/bin/bash

yum update

# Python
yum install nginx redis rabbitmq-server BEDTools
yum install python3
yum install python3-devel libcurl-devel # Needed for pip packages - CyVCF2/pysam

# Postgres 12 (RHEL7 default is quite old)
yum -y install https://download.postgresql.org/pub/repos/yum/reporpms/EL-7-x86_64/pgdg-redhat-repo-latest.noarch.rpm
yum -y install epel-release yum-utils
yum-config-manager --enable pgdg12
yum install postgresql12-server postgresql12 postgresql12-contrib

/usr/pgsql-12/bin/postgresql-12-setup initdb
systemctl enable --now postgresql-12

