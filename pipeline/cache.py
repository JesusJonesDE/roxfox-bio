import json
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any, Optional

from pipeline.config import Settings


class CacheManager:
    def __init__(self, settings: Settings):
        self.settings = settings
        self._manifest: dict = {}
        self._load_manifest()

    # ── Manifest ──────────────────────────────────────────────────────────────

    def _load_manifest(self) -> None:
        path = self.settings.manifest_path
        if path.exists():
            try:
                self._manifest = json.loads(path.read_text())
            except (json.JSONDecodeError, OSError):
                self._manifest = {}
        else:
            self._manifest = {}

    def _save_manifest(self) -> None:
        path = self.settings.manifest_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self._manifest, indent=2))

    # ── Stage completion tracking ──────────────────────────────────────────────

    def is_stage_complete(self, stage: str, target_gene: str) -> bool:
        key = f"{stage}:{target_gene}"
        entry = self._manifest.get(key)
        if not entry:
            return False
        fetched_at = datetime.fromisoformat(entry["fetched_at"])
        age = datetime.now(timezone.utc) - fetched_at
        return age.days < self.settings.cache_max_age_days

    def mark_stage_complete(self, stage: str, target_gene: str) -> None:
        key = f"{stage}:{target_gene}"
        self._manifest[key] = {"fetched_at": datetime.now(timezone.utc).isoformat()}
        self._save_manifest()

    def reset_stage(self, stage: str, target_gene: str) -> None:
        key = f"{stage}:{target_gene}"
        self._manifest.pop(key, None)
        # Cascading resets: resetting fetch also invalidates analyze and report
        cascade = {"fetch": ["analyze", "report"], "analyze": ["report"], "report": []}
        for downstream in cascade.get(stage, []):
            self._manifest.pop(f"{downstream}:{target_gene}", None)
        self._save_manifest()

    def get_stage_date(self, stage: str, target_gene: str) -> Optional[str]:
        key = f"{stage}:{target_gene}"
        entry = self._manifest.get(key)
        return entry["fetched_at"] if entry else None

    # ── Source-level cache (raw API responses) ─────────────────────────────────

    def get_cache_path(self, target_gene: str, source: str) -> Optional[Path]:
        """Return the most recent valid cache file for this target+source, or None."""
        cache_dir = self.settings.cache_dir / target_gene
        if not cache_dir.exists():
            return None
        matches = sorted(cache_dir.glob(f"{source}_*.json"), reverse=True)
        for path in matches:
            if self._is_file_fresh(path):
                return path
        return None

    def _is_file_fresh(self, path: Path) -> bool:
        try:
            stem = path.stem  # e.g. "open_targets_2026-05-26T143000"
            ts_str = stem.rsplit("_", 1)[1].replace("T", " ")
            fetched_at = datetime.fromisoformat(ts_str).replace(tzinfo=timezone.utc)
            age = datetime.now(timezone.utc) - fetched_at
            return age.days < self.settings.cache_max_age_days
        except (IndexError, ValueError):
            return False

    def save(self, target_gene: str, source: str, data: Any, record_count: int) -> Path:
        cache_dir = self.settings.cache_dir / target_gene
        cache_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%S")
        path = cache_dir / f"{source}_{ts}.json"
        path.write_text(json.dumps({"fetched_at": ts, "record_count": record_count, "data": data}, indent=2))
        return path

    def load(self, target_gene: str, source: str) -> Optional[Any]:
        path = self.get_cache_path(target_gene, source)
        if path is None:
            return None
        try:
            payload = json.loads(path.read_text())
            return payload.get("data")
        except (json.JSONDecodeError, OSError):
            return None

    def is_source_fresh(self, target_gene: str, source: str) -> bool:
        return self.get_cache_path(target_gene, source) is not None

    def get_source_meta(self, target_gene: str, source: str) -> Optional[dict]:
        path = self.get_cache_path(target_gene, source)
        if path is None:
            return None
        try:
            payload = json.loads(path.read_text())
            return {"fetched_at": payload.get("fetched_at"), "record_count": payload.get("record_count", 0)}
        except (json.JSONDecodeError, OSError):
            return None
