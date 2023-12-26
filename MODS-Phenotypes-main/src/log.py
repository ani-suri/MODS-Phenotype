import logging
import os
import sys
from IPython.display import display
from pathlib import Path

class JupyterStreamHandler(logging.StreamHandler):
    """
    Custom logging handler sending logs to IPython's display() method.
    """
    def emit(self, record):
        try:
            msg = self.format(record)
            display(msg)
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logger(name, root_folder=Path.home()):
    logger = logging.getLogger(name)
    logger.setLevel(logging.WARNING)  # Set to WARNING to capture both WARNINGS and ERRORS

    # Console handler to print messages to the console
    # console_handler = logging.StreamHandler()
    console_handler = JupyterStreamHandler()
    console_handler.setLevel(logging.WARNING)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
  
    # File handler to save log messages to a file
    script_name = os.path.basename(sys.argv[0])  # Get the name of the script that's running
    log_filename = os.path.splitext(script_name)[0] + ".log"
    # Use the provided root_folder path
    log_filepath = os.path.join(root_folder, log_filename)
    # log_filepath = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), log_filename)

    # File handler to save log messages to a file
    file_handler = logging.FileHandler(log_filepath)
    file_handler.setLevel(logging.WARNING)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger
