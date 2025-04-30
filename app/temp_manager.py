import tempfile
import shutil
from pathlib import Path
from typing import Optional, Generator
from contextlib import contextmanager
import os
from .utils import setup_logging

logger = setup_logging()

class TempFileManager:
    """Manages temporary files for custom network and pathway uploads in PEANUT.
    
    This class handles the creation and cleanup of temporary files used during
    pipeline execution, particularly for custom network and pathway files.
    
    Attributes:
        _temp_dir: Path to the temporary directory
        network_path: Path to the temporary network file
        pathway_file_path: Path to the temporary pathway file
    """
    
    def __init__(self):
        """Initialize the temporary file manager."""
        self._temp_dir: Optional[str] = None
        self.network_path: Optional[Path] = None
        self.pathway_file_path: Optional[Path] = None
        logger.debug("TempFileManager initialized")
    
    @contextmanager
    def create_temp_files(self, network_filename: Optional[str] = None, pathway_filename: Optional[str] = None) -> Generator['TempFileManager', None, None]:
        """Create temporary files for network and pathway data.
        
        Args:
            network_filename: Name for the network file (if needed)
            pathway_filename: Name for the pathway file (if needed)
            
        Yields:
            Self for use in a context manager
            
        Raises:
            OSError: If temporary directory creation fails
            ValueError: If filenames are invalid
        """
        logger.debug("Creating temporary files")
        try:
            # Create a temporary directory with a specific prefix
            self._temp_dir = tempfile.mkdtemp(prefix="peanut_")
            temp_path = Path(self._temp_dir)
            
            # Set appropriate permissions
            os.chmod(self._temp_dir, 0o700)
            
            # Create temporary network file if provided
            if network_filename:
                if not self._is_valid_filename(network_filename):
                    raise ValueError(f"Invalid network filename: {network_filename}")
                self.network_path = temp_path / network_filename
                self.network_path.touch()
                os.chmod(self.network_path, 0o600)
                logger.debug(f"Temporary network path: {self.network_path}")

            # Create temporary pathway file if provided
            if pathway_filename:
                if not self._is_valid_filename(pathway_filename):
                    raise ValueError(f"Invalid pathway filename: {pathway_filename}")
                self.pathway_file_path = temp_path / pathway_filename
                self.pathway_file_path.touch()
                os.chmod(self.pathway_file_path, 0o600)
                logger.debug(f"Temporary pathway path: {self.pathway_file_path}")

            yield self

        except Exception as e:
            logger.error(f"Error creating temporary files: {e}")
            self.cleanup()
            raise
        finally:
            pass

    def _is_valid_filename(self, filename: str) -> bool:
        """Check if a filename is valid.
        
        Args:
            filename: The filename to check
            
        Returns:
            True if the filename is valid, False otherwise
        """
        # Basic validation - you might want to add more checks
        return bool(filename) and not any(c in filename for c in '/\\')

    def cleanup(self) -> None:
        """Clean up temporary files and directory.
        
        This method ensures that all temporary files and the temporary directory
        are properly removed, even if errors occur during cleanup.
        """
        logger.debug("Starting cleanup")
        try:
            if self._temp_dir and Path(self._temp_dir).exists():
                # First, try to remove individual files
                for path in [self.network_path, self.pathway_file_path]:
                    if path and path.exists():
                        try:
                            path.unlink()
                        except Exception as e:
                            logger.warning(f"Failed to remove temporary file {path}: {e}")
                
                # Then remove the directory
                try:
                    shutil.rmtree(self._temp_dir)
                    logger.info("Temporary files cleaned up successfully")
                except Exception as e:
                    logger.error(f"Failed to remove temporary directory {self._temp_dir}: {e}")
        except Exception as e:
            logger.error(f"Error during cleanup: {e}")
        finally:
            # Reset all attributes
            self._temp_dir = None
            self.network_path = None
            self.pathway_file_path = None
