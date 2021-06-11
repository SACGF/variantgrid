-- Need Case Insensitive Text extension, install globally
CREATE EXTENSION IF NOT EXISTS citext;
CREATE DATABASE snpdb2;
CREATE USER snpdb WITH PASSWORD 'snpdb';
GRANT ALL PRIVILEGES ON DATABASE snpdb2 to snpdb;
ALTER USER snpdb CREATEDB;

-- Install extension install on snpdb database
\connect snpdb2;
CREATE EXTENSION IF NOT EXISTS citext;
