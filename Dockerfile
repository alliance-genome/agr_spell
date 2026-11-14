# syntax=docker/dockerfile:1
# spell-next: one image for every MOD instance. Data lives outside the image, in a
# volume mounted at /srv/spell (raw/, builds/, current, logs/; see spell_next.cli).
FROM python:3.13-slim

ARG AG_GRID_VERSION=36.3.0
ARG AG_GRID_SHA256=5fc8736cc0726dcd2e2901681180462fe4d1a519b31a2b4ab7268a97a52a3532

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    SPELL_ROOT=/srv/spell

WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
COPY configs ./configs
RUN pip install --no-cache-dir . \
    && python - <<'EOF'
# Serve AG Grid from the image instead of a CDN; the checksum pins the exact file.
import hashlib, os, sysconfig, urllib.request
version, expected = os.environ['AG_GRID_VERSION'], os.environ['AG_GRID_SHA256']
url = 'https://cdn.jsdelivr.net/npm/ag-grid-community@%s/dist/ag-grid-community.min.js' % version
data = urllib.request.urlopen(url, timeout=120).read()
if hashlib.sha256(data).hexdigest() != expected:
    raise SystemExit('AG Grid checksum mismatch for %s' % url)
target = os.path.join(sysconfig.get_paths()['purelib'], 'spell_next', 'static', 'vendor')
os.makedirs(target, exist_ok=True)
with open(os.path.join(target, 'ag-grid-community.min.js'), 'wb') as fh:
    fh.write(data)
EOF

RUN useradd --create-home --uid 1000 spell
USER spell

EXPOSE 8000
# One-off containers (validate, rebuild; SPELL_ONEOFF=1) run no web app, so they pass trivially.
HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
    CMD [ "$SPELL_ONEOFF" = "1" ] || python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/info', timeout=4)"
CMD ["spell", "serve", "/srv/spell", "--host", "0.0.0.0", "--port", "8000"]
