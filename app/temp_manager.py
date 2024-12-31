import logging
import tempfile
import shutil
from pathlib import Path
from typing import Optional
from contextlib import contextmanager

logger = logging.getLogger(__name__)

class TempFileManager:
    """Manages temporary files for custom network and pathway uploads."""
    
    def __init__(self):
        self._temp_dir: Optional[str] = None
        self.network_path: Optional[Path] = None
        self.pathway_file_path: Optional[Path] = None
        logger.debug("TempFileManager initialized")
    
    @contextmanager
    def create_temp_files(self, network_filename: Optional[str] = None, pathway_filename: Optional[str] = None):
        logger.debug("Creating temporary files")
        try:
            # Create a temporary directory
            self._temp_dir = tempfile.mkdtemp(prefix="peanut_")
            temp_path = Path(self._temp_dir)

            # Create temporary network file if provided
            if network_filename:
                self.network_path = temp_path / network_filename
                logger.debug(f"Temporary network path: {self.network_path}")

            # Create temporary pathway file if provided
            if pathway_filename:
                self.pathway_file_path = temp_path / pathway_filename
                logger.debug(f"Temporary pathway path: {self.pathway_file_path}")

            yield self  # Provide the instance to the caller

        except Exception as e:
            logger.error(f"Error creating temporary files: {e}")
            raise
        finally:
            self.cleanup()

    def cleanup(self):
        """Cleans up temporary files and directory."""
        logger.debug("Starting cleanup")
        try:
            if self._temp_dir and Path(self._temp_dir).exists():
                shutil.rmtree(self._temp_dir)
                logger.info("Temporary files cleaned up successfully")
        except Exception as e:
            logger.error(f"Error cleaning up temporary files: {e}")
        finally:
            self._temp_dir = None
            self.network_path = None
            self.pathway_file_path = None
