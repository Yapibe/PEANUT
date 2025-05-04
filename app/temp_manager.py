import tempfile
import shutil
from pathlib import Path
from typing import Optional, Generator
from contextlib import contextmanager
import os

class TempFileManager:
    """Manages temporary files for custom network and pathway uploads in PEANUT."""
    
    def __init__(self):
        self._temp_dir: Optional[str] = None
        self.network_path: Optional[Path] = None
        self.pathway_file_path: Optional[Path] = None
    
    @contextmanager
    def create_temp_files(self, network_filename: Optional[str] = None, pathway_filename: Optional[str] = None) -> Generator['TempFileManager', None, None]:
        """Create temporary files for network and pathway data."""
        try:
            self._temp_dir = tempfile.mkdtemp(prefix="peanut_")
            temp_path = Path(self._temp_dir)
            os.chmod(self._temp_dir, 0o700)
            
            if network_filename:
                if not self._is_valid_filename(network_filename):
                    raise ValueError(f"Invalid network filename: {network_filename}")
                self.network_path = temp_path / network_filename
                self.network_path.touch()
                os.chmod(self.network_path, 0o600)

            if pathway_filename:
                if not self._is_valid_filename(pathway_filename):
                    raise ValueError(f"Invalid pathway filename: {pathway_filename}")
                self.pathway_file_path = temp_path / pathway_filename
                self.pathway_file_path.touch()
                os.chmod(self.pathway_file_path, 0o600)

            yield self

        except Exception:
            self.cleanup()
            raise

    def _is_valid_filename(self, filename: str) -> bool:
        """Check if a filename is valid."""
        return bool(filename) and not any(c in filename for c in '/\\')

    def cleanup(self) -> None:
        """Clean up temporary files and directory."""
        try:
            if self._temp_dir and Path(self._temp_dir).exists():
                for path in [self.network_path, self.pathway_file_path]:
                    if path and path.exists():
                        path.unlink()
                shutil.rmtree(self._temp_dir)
        finally:
            self._temp_dir = None
            self.network_path = None
            self.pathway_file_path = None
