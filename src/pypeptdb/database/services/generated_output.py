from bson.objectid import ObjectId
from datetime import datetime, timezone
from pathlib import Path

from pypeptdb.database.utils.db_connection import get_db
from pypeptdb.log import get_logger
logger = get_logger(__name__)


db = get_db()
output_collection = db["generated_structures"]


def save_structure_metadata(sequence_biln: str, smiles: str, pdb: str, sdf: str, pdb_filename: str, zip_filename: str, user_id: str=None, pdf_filename: str, note: str=None) -> str:
    """Store the metadata of a generated structure in the database.

    Args:
        sequence_biln (str): BILN sequence
        smiles (str): SMILES representation
        pdb (str): PDB string content
        sdf (str): SDF string content
        pdb_filename (str): PDB filename
        zip_filename (str): ZIP filename
        pdf_filename (str, optional): PDF report filename
        user_id (str, optional): User ID. Defaults to None.
        note (str, optional): Note to describe the output. Defaults to None.
    Returns:
        str: ID of the inserted document
    """
    
    try:
        doc = {
            "sequence_biln": sequence_biln,
            "created_at": datetime.now(timezone.utc),
            "smiles": smiles,
            "pdb": pdb,
            "sdf": sdf,
            "pdf_filename": pdf_filename,
            "pdb_filename": pdb_filename,
            "zip_path": zip_filename,
            "has_pdf": bool(pdf_filename),
            "user_id": user_id,  # will be populated later
            "note": None,  # Placeholder for future use
        }
        result = output_collection.insert_one(doc)
        return str(result.inserted_id)
    except Exception as e:
        logger.error(f"Failed to save structure metadata: {e}")
        raise e


def get_structure_by_id(doc_id: str) -> dict | None:
    """Retrieve a structure document from the database by its ID.

    Args:
        doc_id (str): The ID of the document to retrieve.

    Returns:
        dict | None: The document if found, otherwise None.
    """
    try:
        return output_collection.find_one({"_id": ObjectId(doc_id)})
    except Exception as e:
        logger.error(f"Failed to retrieve structure by ID {doc_id}: {e}")
        return None
    
