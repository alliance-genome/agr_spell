"""Refresh a MOD's reference files (gene names, GO) from their public download URLs.

``mod.toml`` lists the files to fetch, as paths relative to the raw directory::

    [updates]
    "meta/SGD_features.tab" = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"

Each file is downloaded next to its target, checked, and only then swapped in; the file
it replaces is kept as ``<name>.prev``. A download that fails a check never replaces
anything: it must be non-empty, not an HTML error page, a complete gzip stream for
``.gz`` files, and at least MIN_SIZE_RATIO of the current file's size (a truncated or
emptied upstream file would otherwise wipe out gene names or GO annotations).
"""
import gzip
import os
import shutil
import urllib.request
import zlib

MIN_SIZE_RATIO = 0.5
TIMEOUT = 300
USER_AGENT = 'spell-next reference update'


class RefError(Exception):
    pass


def _download(url, dest):
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp, open(dest, 'wb') as out:
        shutil.copyfileobj(resp, out)


def _check(new, current):
    size = os.path.getsize(new)
    if size == 0:
        raise RefError('download is empty')
    if new.endswith('.gz.download'):
        try:
            with gzip.open(new, 'rb') as fh:
                head = fh.read(512)
                while fh.read(1 << 20):
                    pass
        except (OSError, EOFError, zlib.error) as err:
            raise RefError('not a complete gzip file: %s' % err)
    else:
        with open(new, 'rb') as fh:
            head = fh.read(512)
    if head.lstrip().lower().startswith((b'<!doctype html', b'<html')):
        raise RefError('download is an HTML page, not a data file')
    if os.path.exists(current) and size < MIN_SIZE_RATIO * os.path.getsize(current):
        raise RefError('download is %d bytes, under %d%% of the current file (%d bytes)'
                       % (size, MIN_SIZE_RATIO * 100, os.path.getsize(current)))


def _same(a, b):
    if not os.path.exists(b) or os.path.getsize(a) != os.path.getsize(b):
        return False
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        while True:
            x, y = fa.read(1 << 20), fb.read(1 << 20)
            if x != y:
                return False
            if not x:
                return True


def update(raw_dir, updates, log=print):
    """Fetch every file in ``updates`` ({relative path: url}); return (changed, failed) path lists."""
    changed, failed = [], []
    for rel, url in updates.items():
        target = os.path.join(raw_dir, rel)
        tmp = target + '.download'
        try:
            os.makedirs(os.path.dirname(target), exist_ok=True)
            _download(url, tmp)
            _check(tmp, target)
            if _same(tmp, target):
                os.remove(tmp)
                log('  unchanged  %s' % rel)
                continue
            if os.path.exists(target):
                os.replace(target, target + '.prev')
            os.replace(tmp, target)
            changed.append(rel)
            log('  UPDATED    %s (%d bytes; previous kept as %s.prev)' % (rel, os.path.getsize(target), rel))
        except (OSError, RefError) as err:
            if os.path.exists(tmp):
                os.remove(tmp)
            failed.append(rel)
            log('  FAILED     %s: %s (current file kept)' % (rel, err))
    return changed, failed
