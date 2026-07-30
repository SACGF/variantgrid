These service files work for any install location - copy them to /lib/systemd/system as-is.

Set where VariantGrid is installed once, in config/variantgrid.env (copied to
/etc/variantgrid/variantgrid.env), which every service reads:

    VG_INSTALL_DIR="/opt/variantgrid"

Everything else is relative to that - the services cd there before starting, and the
celery/gunicorn .env files run out of .venv/bin (see the "Install system services" wiki page).
