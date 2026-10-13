These service files work for any install location - copy them to /lib/systemd/system as-is.

Set where VariantGrid is installed once, in config/variantgrid.env (copied to
/etc/variantgrid/variantgrid.env), which every service reads:

    VG_INSTALL_DIR="/opt/variantgrid"

Everything else is relative to that - the services cd there before starting, and the
celery/gunicorn .env files run out of .venv/bin (see the "Install system services" wiki page).

The celery units set OOMPolicy=continue: under systemd's default (stop) a kernel OOM kill of one
subprocess (eg somalier) stops the whole unit, SIGTERMing every in-flight task on that queue.
After copying the files: sudo systemctl daemon-reload, then restart the celery units.

redis-server.service.d/variantgrid.conf is a drop-in for the distribution's redis unit (not a unit of
ours). It lowers redis's OOM score so the kernel kills a celery worker instead. The drop-in directory is
named after the unit, which differs by distribution - check with `systemctl list-units 'redis*'`:

    Ubuntu/Debian (redis-server.service):
        sudo mkdir -p /etc/systemd/system/redis-server.service.d
        sudo cp config/systemd/redis-server.service.d/variantgrid.conf /etc/systemd/system/redis-server.service.d/
        sudo systemctl daemon-reload && sudo systemctl restart redis-server

    RHEL (redis.service):
        sudo mkdir -p /etc/systemd/system/redis.service.d
        sudo cp config/systemd/redis-server.service.d/variantgrid.conf /etc/systemd/system/redis.service.d/
        sudo systemctl daemon-reload && sudo systemctl restart redis

Verify with `systemctl show redis -p OOMScoreAdjust` (expect -900).
