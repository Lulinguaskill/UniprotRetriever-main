CREATE TABLE YourTableName (
    uniprot_id VARCHAR(255),
    title VARCHAR(255),
    length VARCHAR(255),
    sequence VARCHAR(MAX),
    e_value FLOAT,
    PRIMARY KEY (uniprot_id, title)
);