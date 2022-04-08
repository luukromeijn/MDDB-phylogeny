import mysql.connector
from mysql.connector import Error

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

        cursor = self.db_connection.cursor()
        try:
            cursor.execute(query)
            if commit:
                self.db_connection.commit()
            result = cursor.fetchall()
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