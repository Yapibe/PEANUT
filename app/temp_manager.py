import logging
import tempfile
import shutil
from pathlib import Path
from typing import Optional
from contextlib import contextmanager

logger = logging.getLogger(__name__)

class TempFileManager:
    """Manages temporary files for custom network uploads."""
    
    def __init__(self):
        self._temp_dir: Optional[str] = None
        self.network_path: Optional[Path] = None
        self.matrix_path: Optional[Path] = None
        logger.debug("TempFileManager initialized")
    
    @contextmanager
    def create_temp_files(self, network_filename: str):
        """
        Creates temporary directory and files for a custom network upload.
        
        Args:
            network_filename (str): Original filename of the uploaded network
        """
        logger.debug(f"Starting create_temp_files for {network_filename}")
        try:
            # Create temporary directory using mkdtemp
            self._temp_dir = tempfile.mkdtemp(prefix="peanut_")
            temp_path = Path(self._temp_dir)
            logger.info(f"Created temporary directory: {temp_path}")
            
            # Create paths for network and matrix files
            self.network_path = temp_path / network_filename
            matrix_name = f"{Path(network_filename).stem}_matrix.npz"
            self.matrix_path = temp_path / matrix_name
            
            logger.debug(f"Temporary network path: {self.network_path}")
            logger.debug(f"Temporary matrix path: {self.matrix_path}")
            
            logger.debug("About to yield TempFileManager instance")
            yield self
            logger.debug("After yield in create_temp_files")
            
        except Exception as e:
            logger.error(f"Error in create_temp_files: {e}", exc_info=True)
            self.cleanup()
            raise
        finally:
            logger.debug("Exiting create_temp_files context")
            
    def cleanup(self):
        """Cleans up temporary files and directory."""
        logger.debug("Starting cleanup")
        try:
            if self._temp_dir:
                logger.debug(f"Temp directory exists: {self._temp_dir}")
                if Path(self._temp_dir).exists():
                    logger.debug(f"About to remove directory: {self._temp_dir}")
                    shutil.rmtree(self._temp_dir)
                    logger.info("Temporary files cleaned up successfully")
                else:
                    logger.warning(f"Temp directory doesn't exist during cleanup: {self._temp_dir}")
            else:
                logger.debug("No temp directory to clean up")
            
            self._temp_dir = None
            self.network_path = None
            self.matrix_path = None
            logger.debug("Cleanup completed")
        except Exception as e:
            logger.error(f"Error during cleanup: {e}", exc_info=True)
            raise

    def __del__(self):
        """Ensures cleanup on object destruction."""
        logger.debug("TempFileManager destructor called")
        self.cleanup() 