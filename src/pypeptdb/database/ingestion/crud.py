from pathlib import Path
from pypeptdb.database.ingestion.data_cleaners import clean_dataframe
from pypeptdb.database.ingestion.data_serializers import serialize_to_pypeptdb_collections, collect_sdf_document, create_mongodb_indexes
from pypeptdb.database.utils.db_connection import get_db

from pypeptdb.utils.data_pypept import load_sdf_data, read_sdf_file

from pypeptdb.log import get_logger
logger = get_logger(__name__)


def add_monomer_to_pypeptdb(sdf_content: str):
    logger.info('Starting ingestion...')

    try:        
        # Collect/Parse SDF data
        logger.info('Step 1: Parsing SDF data into a pypeptdb document...')
        sdf_document = collect_sdf_document(sdf_content=sdf_content)
        residue = sdf_document[0]["symbol"]

        # Insert SDF data document to db
        db = get_db()
        logger.info(f'Step 2: Adding of SDF monomer data to {db.name} as "global_sdf" collection...')
        db["global_sdf"].insert_one(sdf_document[0])

        # Load SDF from db ; required to compute monomers properties
        logger.info(f'Step 3: Loading SDF data from {db.name}...')
        df = load_sdf_data(from_db=True, residues=[residue])
        
        # Transform df (add SMILES, images)
        logger.info('Step 4: Transform dataframe...')
        df_transformed = clean_dataframe(df)

        # Collect global monomers and compute their properties
        logger.info('Step 5: Compute monomers properties...')
        pypeptdb_collections = serialize_to_pypeptdb_collections(df=df_transformed)

        # Insert global monomers collection to db
        logger.info(f'Step 6: Ingestion of global monomers collection to {db.name}...')
        db["global_monomers"].insert_one(pypeptdb_collections["monomers"][0])

        # Insert monomers image collection to db
        logger.info(f'Step 7: Ingestion of monomers image collection to {db.name}...')
        db["monomer_images"].insert_one(pypeptdb_collections["images"][0])

        # Insert monomers properties collection to db
        logger.info(f'Step 8: Ingestion of monomers properties collection to {db.name}...')
        db["global_properties"].insert_one(pypeptdb_collections["properties"][0])
    except Exception as e:
        logger.error(f'Error during ingestion: {e}')
        raise e
    else:
        logger.info('Ingestion completed successfully.')
        return True