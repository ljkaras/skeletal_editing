import pandas as pd
import sqlite3
import sys
import time

time_start = time.time()

if len(sys.argv) != 2:
    print("Usage: csv_to_db.py <filename>")
    sys.exit(1)

filename = sys.argv[1]
print(f"Reading {filename}")

clean_name = filename.split(".")[0]

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(filename)

# Connect to a SQLite database (this will create the database file if it doesn't exist)
conn = sqlite3.connect(f"{clean_name}.db")

# Write the data from the DataFrame to a SQL table
df.to_sql(clean_name, conn, if_exists='replace', index=False)

# Create an index on the 'SMILES' column
cursor = conn.cursor()
cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_smiles ON {clean_name}(SMILES)")
conn.commit()

# Close the connection to the database
conn.close()

print(f"Database {clean_name}.db created with index on 'SMILES' column.")

time_end = time.time()
print(f"Processing time: {time_end - time_start:.2f} s")
