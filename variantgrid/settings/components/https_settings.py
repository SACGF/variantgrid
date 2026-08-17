"""
Settings that only make sense when a deployment is served over HTTPS.

Imported by the env settings files rather than default_settings, as intranet deployments served over
plain HTTP need these left off - a Secure cookie is never stored by the browser, so marking them
Secure on an HTTP site stops anyone logging in.
"""

SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
