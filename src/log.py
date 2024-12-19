import logging
import sys

def configure_logger():
    """
    Configure the shared logging for the entire application.
    This function is called once and sets the logging configuration globally.
    """
    logging.basicConfig(
        level=logging.INFO,  # Set the global log level
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),  # Print logs to stdout
        ]
    )

def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance with the shared configuration.
    
    Args:
        name (str): The name of the logger (typically use __name__).
    
    Returns:
        logging.Logger: Configured logger instance.
    """
    return logging.getLogger(name)

# Call configure_logger once when this file is imported
configure_logger()