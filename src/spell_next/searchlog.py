"""Monthly search logs, like SPELL's SPELL_NIH_Log_YYYY_MM.log (usage reporting).

One tab-separated line per search: date, time, kind (single | multi | listing | embed),
query genes (systematic names, comma-separated) and dataset tags if any. No client
addresses or other personal data are recorded.
"""
import threading
from datetime import datetime
from pathlib import Path

_lock = threading.Lock()


class SearchLog:
    def __init__(self, directory):
        self.directory = Path(directory) if directory else None

    def record(self, kind, genes, tags=()):
        if self.directory is None:
            return
        now = datetime.now()
        line = '\t'.join([now.strftime('%Y-%m-%d'), now.strftime('%H:%M:%S'), kind,
                          ','.join(genes), '|'.join(tags)]) + '\n'
        path = self.directory / ('searches_%s.tsv' % now.strftime('%Y_%m'))
        with _lock:
            self.directory.mkdir(parents=True, exist_ok=True)
            with open(path, 'a') as fh:
                fh.write(line)
