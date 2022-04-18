import mysql.connector
from mysql.connector import Error
from alfpy.utils import seqrecords
from alfpy import wmetric
from alfpy.utils.data import subsmat
from Chunks import Chunk
import numpy as np

class PairwiseDistance:
    '''Calculates pairwise distances, stores in and retrieves from database'''

    def __init__(self, host_name: str="localhost", user_name: str="mddb-phylogeny", user_password: str="mddb", db_name: str="UniteDistances"):
        self.db_link = ConnectDatabase(host_name, user_name, user_password, db_name)

    
    def calc_from_fasta(self, path: str):
        '''Calculates w-metric pairwise distances and stores in SQL database'''

        matrix = subsmat.get('blosum62')
        if self.db_link.table_exists("distances"): 
            cont = input("This will overwrite the current database. Continue? Y/N ")
            if cont.lower().strip() != "y":
                print("Aborted.")
                return

        # Calculate distances
        self.db_link.create_table("distances", "seq_1 INT, seq_2 INT, distance FLOAT(24), INDEX(seq_1), INDEX(seq_2)")
        fh = open(path)
        seq_records = seqrecords.read_fasta(fh)
        fh.close()
        distances = wmetric.Distance(seq_records, matrix)
        for seq_1 in range(seq_records.count): # NOTE Takes 3728 secondes for 500, should be seq_records.count
            fill_table_query = "INSERT INTO distances VALUES "
            for seq_2 in range(seq_1, seq_records.count):
                fill_table_query += "(" + str(seq_1) + ", " + str(seq_2) + ", " + str(distances.pairwise_distance(seq_1, seq_2)) + "),\n"
            fill_table_query = fill_table_query[:-2] + ";"
            self.db_link.query(fill_table_query, commit=True)
    

    def distance(self, seq_1: int, seq_2: int) -> float:
        '''Returns pairwise distance between two sequences'''
        
        sequences = [seq_1, seq_2]
        sequences.sort()
        query = "SELECT distance FROM distances WHERE seq_1=" + str(sequences[0]) + " AND seq_2=" + str(sequences[1]) + " LIMIT 1;"
        
        return self.db_link.query(query)[0][0]


    def determine_outgroup(self, chunk: Chunk, n: int=10) -> 'list[int]':
        '''Determines outgroup of size n to chunk based on distance to ingroup'''

        # Storing the ingroup as a table
        self.db_link.create_table("ingroup", "seq int, INDEX(seq)")
        ingroup_creation = "INSERT INTO ingroup VALUES "
        for seq_id in chunk.ingroup:
            ingroup_creation += "(" + str(seq_id) + "),"
        self.db_link.query(ingroup_creation[:-1] + ";", commit=True)

        # Retrieving sequence ids
        indices = self.db_link.query("SELECT DISTINCT seq_1 FROM distances")
        indices = [i[0] for i in indices]

        dist_to_outgroup = []
        for i in indices: 
            if i in chunk.ingroup: # Should not be in ingroup
                dist_to_outgroup.append(np.NaN) # Append NaN to keep index structure for argsort
                continue
            query = """
            SELECT AVG(distance) 
            FROM distances 
            WHERE
            (seq_2 in (select seq from ingroup where seq=seq_2) and seq_1=""" + str(i) + """) OR
            (seq_1 in (select seq from ingroup where seq=seq_1) and seq_2=""" + str(i) + """);"""
            dist_to_outgroup.append(self.db_link.query(query)[0][0])

        sorted = np.argsort(dist_to_outgroup)[:n]
        # print([dist_to_outgroup[i] for i in sorted]) # Uncomment to print distances
        outgroup = [seq_index for seq_index in sorted]
        return outgroup


class ConnectDatabase:
    '''Connects and queries a mysql database'''

    def __init__(self, host_name: str, user_name: str, user_password: str, db_name: str):

        self.sql_connection = self.connect(host_name, user_name, user_password)
        if not self.db_exists(db_name):
            self.create_db(db_name)
        self.db_name = db_name
        self.db_connection = self.connect(host_name, user_name, user_password, db_name)

    
    def query(self, query: str, commit: bool=False):
        '''Performs a query on the connected database'''

        self.db_connection.reconnect()
        cursor = self.db_connection.cursor()
        try:
            cursor.execute(query)
            if commit:
                self.db_connection.commit()
            result = cursor.fetchall()
            cursor.close()
            return result
        except Error as err:
            print(f"Error: '{err}'")

    
    def connect(self, host_name: str, user_name: str, user_password: str, db_name: str=None):
        '''Connects to mysql (specific database if specified)'''

        try:
            connection = mysql.connector.connect(
                host=host_name,
                user=user_name,
                passwd=user_password,
                database=db_name
            )
        except Error as err:
            raise RuntimeError(f"Error: '{err}'")

        return connection


    def db_exists(self, db_name: str):
        '''Checks if database exists'''

        cursor = self.sql_connection.cursor()
        query = """
        SELECT COUNT(SCHEMA_NAME)
        FROM INFORMATION_SCHEMA.SCHEMATA
        WHERE SCHEMA_NAME = '""" + db_name + """'
        """
        try:
            cursor.execute(query)
            result = cursor.fetchall()
            return bool(result[0][0])
        except Error as err:
            print(f"Error: '{err}'")
    

    def create_db(self, db_name: str):
        '''Creates an empty database'''

        self.high_level_query("DROP DATABASE IF EXISTS " + db_name)
        self.high_level_query("CREATE DATABASE " + db_name)


    def table_exists(self, table_name: str):
        cursor = self.sql_connection.cursor()
        query = """
        SELECT COUNT(TABLE_NAME)
        FROM INFORMATION_SCHEMA.TABLES
        WHERE TABLE_NAME = '""" + table_name + """'
        """
        try:
            cursor.execute(query)
            result = cursor.fetchall()
            return bool(result[0][0])
        except Error as err:
            print(f"Error: '{err}'")


    def create_table(self, table_name: str, columns: str):
        '''Creates an empty table'''

        self.query("DROP TABLE IF EXISTS " + table_name, commit=True)
        creation_query = "CREATE TABLE " + table_name + " (" + columns + ");"
        self.query(creation_query, commit=True)


    def high_level_query(self, query: str):
        '''Performs a query on sql level'''

        cursor = self.sql_connection.cursor()
        try:
            cursor.execute(query)
        except Error as err:
            print(f"Error: '{err}'")