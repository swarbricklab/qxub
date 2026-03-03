"""
Platform definition loader with URL support.

Loads platform definitions from various sources:
- file:// - Local file paths
- https:// - Remote HTTP/HTTPS URLs with caching
- Relative paths - Search standard locations

Caching is applied only to remote (https://) URLs to avoid
repeated network requests.
"""

import hashlib
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import urlparse

import requests
import yaml

logger = logging.getLogger(__name__)


class PlatformLoadError(Exception):
    """Raised when platform definition cannot be loaded."""

    pass


class PlatformLoader:
    """Load platform definitions from various sources with caching."""

    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl_hours: int = 24):
        """
        Initialize platform loader.

        Args:
            cache_dir: Directory for cached platform definitions
            cache_ttl_hours: Time-to-live for cached definitions in hours
        """
        self.cache_dir = cache_dir or (Path.home() / ".cache" / "qxub" / "platforms")
        self.cache_ttl = timedelta(hours=cache_ttl_hours)

    def load_platform_definition(
        self, definition_url: str, platform_search_paths: Optional[list[Path]] = None
    ) -> dict:
        """
        Load platform definition from URL or path.

        Args:
            definition_url: URL or path to platform definition
            platform_search_paths: Optional list of paths to search for relative paths

        Returns:
            Platform definition dictionary

        Raises:
            PlatformLoadError: If definition cannot be loaded
        """
        parsed = urlparse(definition_url)

        if parsed.scheme == "https" or parsed.scheme == "http":
            return self._load_from_http(definition_url)
        elif parsed.scheme == "file":
            return self._load_from_file(Path(parsed.path))
        elif parsed.scheme == "":
            # No scheme - treat as relative path and search
            return self._load_from_search_paths(definition_url, platform_search_paths)
        else:
            raise PlatformLoadError(
                f"Unsupported URL scheme: {parsed.scheme} in {definition_url}"
            )

    def _load_from_http(self, url: str) -> dict:
        """
        Load platform definition from HTTP(S) URL with caching.

        Args:
            url: HTTP(S) URL to platform definition

        Returns:
            Platform definition dictionary

        Raises:
            PlatformLoadError: If download or parsing fails
        """
        # Check cache first
        cached_def, cached_meta = self._load_from_cache(url)
        if cached_def and self._is_cache_valid(cached_meta):
            logger.debug(f"Using cached platform definition from {url}")
            return cached_def

        # Fetch from URL
        logger.info(f"Fetching platform definition from {url}")
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            # Parse YAML
            definition = yaml.safe_load(response.text)
            if not isinstance(definition, dict):
                raise PlatformLoadError(
                    f"Platform definition must be a YAML dictionary: {url}"
                )

            # Cache the definition
            self._save_to_cache(url, definition, response.headers.get("ETag"))

            return definition

        except requests.RequestException as e:
            # If fetch fails but we have a cached version, use it with a warning
            if cached_def:
                logger.warning(f"Failed to fetch {url}, using stale cache: {e}")
                return cached_def
            raise PlatformLoadError(
                f"Failed to fetch platform definition from {url}: {e}"
            )
        except yaml.YAMLError as e:
            raise PlatformLoadError(f"Invalid YAML in platform definition {url}: {e}")

    def _load_from_file(self, file_path: Path) -> dict:
        """
        Load platform definition from local file.

        Args:
            file_path: Path to platform definition file

        Returns:
            Platform definition dictionary

        Raises:
            PlatformLoadError: If file cannot be read or parsed
        """
        try:
            expanded_path = file_path.expanduser().resolve()
            logger.debug(f"Loading platform definition from {expanded_path}")

            if not expanded_path.exists():
                raise PlatformLoadError(
                    f"Platform definition file not found: {expanded_path}"
                )

            with open(expanded_path, "r", encoding="utf-8") as f:
                definition = yaml.safe_load(f)

            if not isinstance(definition, dict):
                raise PlatformLoadError(
                    f"Platform definition must be a YAML dictionary: {expanded_path}"
                )

            return definition

        except yaml.YAMLError as e:
            raise PlatformLoadError(
                f"Invalid YAML in platform definition {file_path}: {e}"
            )
        except IOError as e:
            raise PlatformLoadError(
                f"Cannot read platform definition file {file_path}: {e}"
            )

    def _load_from_search_paths(
        self, filename: str, search_paths: Optional[list[Path]] = None
    ) -> dict:
        """
        Load platform definition by searching standard locations.

        Args:
            filename: Filename to search for
            search_paths: Optional list of paths to search

        Returns:
            Platform definition dictionary

        Raises:
            PlatformLoadError: If file not found in any search path
        """
        if search_paths is None:
            search_paths = self._get_default_search_paths()

        logger.debug(f"Searching for platform definition: {filename}")
        for search_path in search_paths:
            candidate = search_path / filename
            logger.debug(f"  Checking {candidate}")
            if candidate.exists():
                logger.debug(f"  Found at {candidate}")
                return self._load_from_file(candidate)

        # Not found in any search path
        search_paths_str = ", ".join(str(p) for p in search_paths)
        raise PlatformLoadError(
            f"Platform definition '{filename}' not found in search paths: {search_paths_str}"
        )

    def _get_default_search_paths(self) -> list[Path]:
        """Get default platform definition search paths."""
        return [
            Path("/etc/xdg/qxub/platforms"),
            Path.home() / ".config" / "qxub" / "platforms",
            Path.home() / ".qxub" / "platforms",
        ]

    def _get_cache_path(self, url: str) -> Tuple[Path, Path]:
        """
        Get cache file paths for a URL.

        Args:
            url: URL to generate cache paths for

        Returns:
            Tuple of (definition_path, metadata_path)
        """
        # Create a hash of the URL for the cache filename
        url_hash = hashlib.sha256(url.encode()).hexdigest()[:16]
        # Also use a safe version of the URL for human readability
        safe_url = url.replace("://", "_").replace("/", "_").replace(":", "_")
        if len(safe_url) > 100:
            safe_url = safe_url[:100]

        cache_file = self.cache_dir / f"{safe_url}_{url_hash}.yaml"
        cache_meta = self.cache_dir / f"{safe_url}_{url_hash}.meta.yaml"

        return cache_file, cache_meta

    def _load_from_cache(self, url: str) -> Tuple[Optional[dict], Optional[dict]]:
        """
        Load platform definition and metadata from cache.

        Args:
            url: URL of the platform definition

        Returns:
            Tuple of (definition, metadata) or (None, None) if not cached
        """
        cache_file, cache_meta = self._get_cache_path(url)

        if not cache_file.exists() or not cache_meta.exists():
            return None, None

        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                definition = yaml.safe_load(f)

            with open(cache_meta, "r", encoding="utf-8") as f:
                metadata = yaml.safe_load(f)

            return definition, metadata

        except Exception as e:
            logger.warning(f"Failed to load from cache: {e}")
            return None, None

    def _save_to_cache(self, url: str, definition: dict, etag: Optional[str] = None):
        """
        Save platform definition to cache.

        Args:
            url: URL of the platform definition
            definition: Platform definition dictionary
            etag: Optional ETag from HTTP response
        """
        cache_file, cache_meta = self._get_cache_path(url)

        # Create cache directory if it doesn't exist
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Save definition
            with open(cache_file, "w", encoding="utf-8") as f:
                yaml.dump(definition, f, default_flow_style=False, sort_keys=False)

            # Save metadata
            metadata = {
                "source": url,
                "fetched": datetime.now().isoformat(),
                "etag": etag,
            }
            with open(cache_meta, "w", encoding="utf-8") as f:
                yaml.dump(metadata, f, default_flow_style=False)

            logger.debug(f"Cached platform definition to {cache_file}")

        except Exception as e:
            logger.warning(f"Failed to save to cache: {e}")

    def _is_cache_valid(self, metadata: Optional[dict]) -> bool:
        """
        Check if cached definition is still valid based on TTL.

        Args:
            metadata: Cache metadata dictionary

        Returns:
            True if cache is valid, False otherwise
        """
        if not metadata or "fetched" not in metadata:
            return False

        try:
            fetched = datetime.fromisoformat(metadata["fetched"])
            age = datetime.now() - fetched
            is_valid = age < self.cache_ttl

            if not is_valid:
                logger.debug(f"Cache expired (age: {age}, TTL: {self.cache_ttl})")

            return is_valid

        except (ValueError, TypeError) as e:
            logger.warning(f"Invalid cache metadata: {e}")
            return False

    def clear_cache(self, url: Optional[str] = None):
        """
        Clear cached platform definitions.

        Args:
            url: Optional specific URL to clear. If None, clears all cache.
        """
        if url:
            # Clear specific URL
            cache_file, cache_meta = self._get_cache_path(url)
            if cache_file.exists():
                cache_file.unlink()
                logger.info(f"Cleared cache for {url}")
            if cache_meta.exists():
                cache_meta.unlink()
        else:
            # Clear all cache
            if self.cache_dir.exists():
                for cache_file in self.cache_dir.glob("*.yaml"):
                    cache_file.unlink()
                logger.info(f"Cleared all platform cache from {self.cache_dir}")

    def refresh_cache(self, url: str) -> dict:
        """
        Force refresh of cached platform definition.

        Args:
            url: URL to refresh

        Returns:
            Refreshed platform definition

        Raises:
            PlatformLoadError: If refresh fails
        """
        # Clear cache for this URL
        self.clear_cache(url)

        # Re-fetch
        return self._load_from_http(url)


# Global platform loader instance
platform_loader = PlatformLoader()
