import mysql.connector
from mysql.connector import Error
from alfpy.utils import seqrecords
from alfpy import word_pattern
from alfpy import word_vector
from alfpy import word_distance

def create_server_connection(host_name, user_name, user_password):
    connection = None
    try:
        connection = mysql.connector.connect(
            host=host_name,
            user=user_name,
            passwd=user_password
        )
        print("MySQL Database connection successful")
    except Error as err:
        print(f"Error: '{err}'")

    return connection

connection = create_server_connection("localhost", "root", "root")

def create_database(connection, query):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        print("Database created successfully")
    except Error as err:
        print(f"Error: '{err}'")

drop_database_query = """
DROP DATABASE IF EXISTS distances
"""
create_database(connection, drop_database_query)
create_database_query = """
CREATE DATABASE distances"""
create_database(connection, create_database_query)

def create_db_connection(host_name, user_name, user_password, db_name):
    connection = None
    try:
        connection = mysql.connector.connect(
            host=host_name,
            user=user_name,
            passwd=user_password,
            database=db_name
        )
        print("MySQL Database connection successful")
    except Error as err:
        print(f"Error: '{err}'")

    return connection

def execute_query(connection, query):
    cursor = connection.cursor()
    try:
        cursor.execute(query)
        connection.commit()
        print("Query successful")
    except Error as err:
        print(f"Error: '{err}'")

def execute_multi_query(connection, sql, val):
    cursor = connection.cursor()
    try:
        cursor.executemany(sql, val)
        connection.commit()
        print("Query successful")
    except Error as err:
        print(f"Error: '{err}'")

drop_distance_table = """
DROP TABLE IF EXISTS distances
 """

create_distance_table = """
CREATE TABLE pw_distances (
  seq_1 INT, seq_2 INT,
  pw_distance FLOAT(24)
  );
 """

connection = create_db_connection("localhost", "root", "root", "distances") # Connect to the Database
execute_query(connection, drop_distance_table) # Execute our defined query
execute_query(connection, create_distance_table) # Execute our defined query

fh = open('data/sh_qiime_release_10.05.2021/sh_refs_qiime_ver8_97_10.05.2021.fasta') #NOTE Linux specific
seq_records = seqrecords.read_fasta(fh)
fh.close()

pattern = word_pattern.create(seq_records.seq_list, word_size=3)
counts = word_vector.Counts(seq_records.length_list, pattern)
distances = word_distance.Distance(counts, 'google')


for seq_1 in range(seq_records.count):
    fill_table_query = "INSERT INTO pw_distances VALUES "
    for seq_2 in range(seq_1, seq_records.count):
        fill_table_query += "(" + str(seq_1) + ", " + str(seq_2) + ", " + str(distances.pairwise_distance(seq_1, seq_2)) + "),\n"
    fill_table_query = fill_table_query[:-2] + ";"
    execute_query(connection, fill_table_query)