from mongodb.ingestion.data_cleaners import clean_dataframe
from mongodb.ingestion.data_serializers import serialize_to_mongodb_collections
from mongodb.utils.db_connection import get_db


from utils.data_pypept import load_sdf_file


# from PIL import Image
# import io

# img_buffer = io.BytesIO(df["image_binary"]["C"])
# img = Image.open(img_buffer)
# img.show()  # Display the image


def ingest():
    # 1. Load source file
    df = load_sdf_file()
    
    # 2. Clean and transform data
    df_clean = clean_dataframe(df)
    
    # 3. Serialize for MongoDB
    mongodb_collections = serialize_to_mongodb_collections(df_clean)
    
    # 4. Insert into MongoDB
    db = get_db()
    monomers_collection = db["global_monomers"]
    properties_collection = db["global_properties"]

    monomers_collection.insert_many(mongodb_collections["monomers"])
    properties_collection.insert_many(mongodb_collections["properties"])



