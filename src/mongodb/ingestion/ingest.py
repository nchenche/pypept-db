from pathlib import Path
from mongodb.ingestion.data_cleaners import clean_dataframe
from mongodb.ingestion.data_serializers import serialize_to_mongodb_collections
from mongodb.utils.db_connection import get_db

from utils.data_pypept import load_sdf_data

from log import get_logger
logger = get_logger(__name__)


# img_buffer = io.BytesIO(df["image_binary"]["C"])
# img = Image.open(img_buffer)
# img.show()  # Display the image


def ingest_data_to_mongodb(source: str|Path):
    logger.info('Starting ingestion...')

    try:
        # 1. Load source file
        df = load_sdf_data(from_file=source)
        logger.info('Step 1 completed: SDF monomers data successfully loaded into dataframe.')
        
        # 2. Clean and transform data
        df_clean = clean_dataframe(df)
        logger.info('Step 2 completed: Dataframe cleaned and transformed.')

        
        # 3. Serialize for MongoDB
        mongodb_collections = serialize_to_mongodb_collections(df_clean)
        logger.info('Step 3 completed: Data collections ready to be ingested.')
        
        # 4. Access and insert data to the db
        db = get_db()
        logger.info(f'Step 4.1 completed: Successful access to the database {db.name}.')

        monomers_collection = db["global_monomers"]
        properties_collection = db["global_properties"]
        sdf_collections = db["global_sdf"]

        monomers_collection.drop()
        monomers_collection.insert_many(mongodb_collections["monomers"])
        logger.info(f'Step 4.2 completed: monomers data successfully added as "global_monomers" collection.')

        properties_collection.drop()
        properties_collection.insert_many(mongodb_collections["properties"])
        logger.info(f'Step 4.3 completed: properties data successfully added as "global_properties" collection.')
        
        sdf_collections.drop()
        sdf_collections.insert_many(mongodb_collections["sdf"])
        logger.info(f'Step 4.4 completed: sdf data successfully added as "global_sdf" collection.')

        logger.info('Ingestion completed.')
    except Exception as e:
        logger.error(f"Ingestion failed: {e}")



