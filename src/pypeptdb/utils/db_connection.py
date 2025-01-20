from pymongo import MongoClient


def get_db():
    """Connect to the pypeptdb instance using URL from .env"""
    import os
    mongo_url = os.getenv("MONGO_URL", "mongodb://localhost:27017")
    client = MongoClient(mongo_url)
    return client["pypeptdb"]