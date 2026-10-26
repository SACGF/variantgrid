CREATE DATABASE snpdb;
CREATE USER snpdb WITH PASSWORD 'snpdb';
GRANT ALL PRIVILEGES ON DATABASE snpdb to snpdb;
ALTER USER snpdb CREATEDB;

\connect snpdb;
GRANT ALL ON SCHEMA public TO snpdb;  -- Needed for Postgres 15+
