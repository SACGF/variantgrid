-- Need Case Insensitive Text extension, install globally
CREATE EXTENSION IF NOT EXISTS citext;
CREATE DATABASE snpdb;
CREATE USER snpdb WITH PASSWORD 'snpdb';
GRANT ALL PRIVILEGES ON DATABASE snpdb to snpdb;
ALTER USER snpdb CREATEDB;

-- Install extension install on snpdb database
\connect snpdb;
CREATE EXTENSION IF NOT EXISTS citext;
