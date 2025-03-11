from io import BytesIO, StringIO
from pathlib import Path
from datetime import datetime, timezone
import pytest

from pypeptdb.database.utils.db_connection import get_db
from pypeptdb.database.ingestion.data_serializers import (
    collect_sdf_document,
    serialize_to_pypeptdb_collections
)
from pypeptdb.utils.data_pypept import load_sdf_data, get_combined_sdf, get_monomers_df

from pypeptdb.tests.db.lib import _clean_dataframe, _serialize_to_pypeptdb_collections

from rdkit import Chem



SDF_MONOMER = """
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.8662    0.2994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0002   -0.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8658    0.3006    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7322   -0.1986    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.7318   -0.2014    0.0000 R#  0  0  0  0  0  1  0  0  0  0  0  0
  5  1  1  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  RGP 1   5   1
V    5 *
M  END

>  <m_name>
Mycap

>  <symbol>
mc

>  <m_abbr>
mc

>  <m_type>
cap

>  <m_subtype>
natural

>  <m_Rgroups>
None,H,None,None

>  <m_RgroupIdx>
None,4,None,None

>  <m_attachmentPointIdx>
None,0,None,None

>  <natAnalog>
X

>  <pdbName>
MCX

$$$$
"""


# Fixture for using a test database with mongomock
@pytest.fixture
def mock_db():
    db = get_db(mock=True)
    return db


def test_db_access(mock_db):
    assert mock_db.name == "test_db"


def test_collect_sdf():
    KEYS_IN_DOC = ['_id', 'symbol', 'sdf', 'created_at']
    sdf = collect_sdf_document(SDF_MONOMER)
    assert len(sdf) == 1
    assert all([ key in sdf[0].keys() for key in KEYS_IN_DOC ])
    assert sdf[0]["symbol"] == "mc"


def test_insert_sdf_to_db(mock_db):
    sdf = collect_sdf_document(SDF_MONOMER)[0]
    mock_db["global_sdf"].drop()
    mock_db["global_sdf"].insert_one(sdf)
    assert mock_db["global_sdf"].count_documents({}) == 1
    assert mock_db["global_sdf"].find_one()["symbol"] == "mc"


def test_combined_sdf(mock_db):
    mock_db["global_sdf"].drop()
    mock_db["global_sdf"].insert_one(collect_sdf_document(SDF_MONOMER)[0])
    combined_sdf = get_combined_sdf(symbols=['mc'], mock=True)
    assert combined_sdf == SDF_MONOMER


def test_load_sdf_from_db():
    df = load_sdf_data(from_db=True, mock=True)

    expected_columns = [
        'm_name', 'm_abbr', 'm_type', 'm_subtype', 
        'm_Rgroups', 'm_RgroupIdx', 'm_attachmentPointIdx', 
        'natAnalog', 'pdbName', 'ID', 'm_romol'
    ]
    assert sorted(list(df.columns)) == sorted(expected_columns)
    # Verify that the index is as expected
    assert 'mc' in df.index

    # Retrieve the row and perform content checks
    row = df.loc['mc']
    assert row['m_name'] == 'Mycap'
    assert row['m_abbr'] == 'mc'
    assert row['m_type'] == 'cap'
    assert row['m_subtype'] == 'natural'

    # Ensure m_romol is an RDKit molecule
    assert isinstance(row['m_romol'], Chem.rdchem.Mol)


def test_clean_dataframe():
    df = load_sdf_data(from_db=True, mock=True)
    df_transformed = _clean_dataframe(df)

    assert 'smiles' in df_transformed.columns
    assert 'canonic_smiles' in df_transformed.columns
    assert 'image_binary' in df_transformed.columns
    assert 'symbol' in df_transformed.columns
    assert '_id' in df_transformed.columns
    assert 'm_romol' not in df_transformed.columns


def test_serialize_to_pypeptdb_collections():
    df = load_sdf_data(from_db=True, mock=True)
    df_transformed = _clean_dataframe(df)

    pypeptdb_collections = _serialize_to_pypeptdb_collections(df_transformed)
    assert 'monomers' in pypeptdb_collections.keys()
    assert pypeptdb_collections['monomers'][0]['symbol'] == 'mc'

    assert 'properties' in pypeptdb_collections.keys()
    assert pypeptdb_collections['properties'][0]['symbol'] == 'mc'

    assert 'images' in pypeptdb_collections.keys()
    assert pypeptdb_collections['images'][0]['symbol'] == 'mc'
    assert pypeptdb_collections['images'][0]['image_binary'] is not None


def test_insert_monomer_to_global_monomers(mock_db):
    df = load_sdf_data(from_db=True, mock=True)
    df_transformed = _clean_dataframe(df)
    pypeptdb_collections = _serialize_to_pypeptdb_collections(df_transformed)

    mock_db["global_monomers"].drop()
    mock_db["global_monomers"].insert_one(pypeptdb_collections["monomers"][0])
    assert mock_db["global_monomers"].count_documents({}) == 1
    assert mock_db["global_monomers"].find_one()["symbol"] == "mc"


def test_insert_images_to_monomer_images(mock_db):
    df = load_sdf_data(from_db=True, mock=True)
    df_transformed = _clean_dataframe(df)
    pypeptdb_collections = _serialize_to_pypeptdb_collections(df_transformed)

    mock_db["monomer_images"].drop()
    mock_db["monomer_images"].insert_one(pypeptdb_collections["images"][0])
    assert mock_db["monomer_images"].count_documents({}) == 1
    assert mock_db["monomer_images"].find_one()["symbol"] == "mc"


def test_insert_properties_to_global_properties(mock_db):
    df = load_sdf_data(from_db=True, mock=True)
    df_transformed = _clean_dataframe(df)
    pypeptdb_collections = _serialize_to_pypeptdb_collections(df_transformed)

    mock_db["global_properties"].drop()
    mock_db["global_properties"].insert_one(pypeptdb_collections["properties"][0])
    assert mock_db["global_properties"].count_documents({}) == 1
    assert mock_db["global_properties"].find_one()["symbol"] == "mc"


