from pymongo import MongoClient


def get_db(mock: bool = False):
    """Connect to the pypeptdb instance using URL from .env"""
    import os
    mongo_url = os.getenv("MONGO_URL", "mongodb://localhost:27017")
    client = MongoClient(mongo_url)
    db_name = "pypeptdb" if mock == False else "test_db"
    return client[db_name]