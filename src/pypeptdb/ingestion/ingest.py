from pathlib import Path
from pypeptdb.ingestion.data_cleaners import clean_dataframe
from pypeptdb.ingestion.data_serializers import serialize_to_pypeptdb_collections, collect_sdf_from_file, create_mongodb_indexes
from pypeptdb.utils.db_connection import get_db

from utils.data_pypept import load_sdf_data

from log import get_logger
logger = get_logger(__name__)


# img_buffer = io.BytesIO(df["image_binary"]["C"])
# img = Image.open(img_buffer)
# img.show()  # Display the image


def ingest_data_to_pypeptdb(source: str|Path):
    logger.info('Starting ingestion...')

    try:
        # 1. Load source file
        logger.info('Step 1: Loading SDF monomers data into dataframe...')
        df = load_sdf_data(from_file=source)
        
        # 2. Collect/Parse SDF data
        logger.info('Step 2: Parsing SDF monomers data...')
        sdf = collect_sdf_from_file(path=source)
        
        # 3. Insert SDF data collection to db
        db = get_db()
        logger.info(f'Step 3: Ingestion of SDF monomers data to {db.name} as "global_sdf" collection...')
        db["global_sdf"].drop()
        db["global_sdf"].insert_many(sdf)

        # 4. Load SDF from db
        logger.info(f'Step 4: Loading SDF data collection from {db.name}...')
        df = load_sdf_data(from_db=True)
        
        # 5. Transform df (add SMILES, images)
        logger.info('Step 5: Transform dataframe...')
        df_transformed = clean_dataframe(df)

        # 6. Collect global monomers and compute their properties
        logger.info('Step 6: Compute monomers properties...')
        pypeptdb_collections = serialize_to_pypeptdb_collections(df=df_transformed)

        # 7. Insert global monomers collection to db
        logger.info(f'Step 7: Ingestion of global monomers collection to {db.name}...')
        db["global_monomers"].drop()
        db["global_monomers"].insert_many(pypeptdb_collections["monomers"])

        # 8. Insert monomers properties collection to db
        logger.info(f'Step 8: Ingestion of monomers properties collection to {db.name}...')
        db["global_properties"].drop()
        db["global_properties"].insert_many(pypeptdb_collections["properties"])

        # 9. Create useful indexes
        logger.info(f'Step 9: Creation of mongodb indexes...')
        create_mongodb_indexes(db=db)


        logger.info('Ingestion completed.')
    except Exception as e:
        logger.error(f"Ingestion failed: {e}")



