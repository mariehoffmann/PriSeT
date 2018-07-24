SET FOREIGN_KEY_CHECKS = 0;
DROP TABLE IF EXISTS `nodes`;
DROP TABLE IF EXISTS `division`;
DROP TABLE IF EXISTS `gencode`;
DROP TABLE IF EXISTS `names`;
DROP TABLE IF EXISTS `delnodes`;
DROP TABLE IF EXISTS `merged`;
DROP TABLE IF EXISTS `citations`;
DROP TABLE IF EXISTS `typeoftype`;
DROP TABLE IF EXISTS `host`;
DROP TABLE IF EXISTS `type_material`;
DROP TABLE IF EXISTS `rankedlineage`;
DROP TABLE IF EXISTS `fullnamelineage`;
SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE `nodes` (
    `tax_id` INTEGER NOT NULL,
    `parent_tax_id` INTEGER NOT NULL,
    `rank` VARCHAR NOT NULL,
    `embl_code` VARCHAR NOT NULL,
    `division_id` INTEGER NOT NULL,
    `inherited_div_flag` BOOL NOT NULL,
    `genetic_code_id` INTEGER NOT NULL,
    `inherited_GC_flag` BOOL NOT NULL,
    `mitochondrial_genetic_code_id` INTEGER NOT NULL,
    `inherited_MGC_flag` BOOL NOT NULL,
    `GenBank_hidden_flag` BOOL NOT NULL,
    `hidden_subtree_root_flag` BOOL NOT NULL,
    `comments` TEXT NOT NULL,
    `plastid_genetic_code_id` INTEGER NOT NULL,
    `inherited_PGC_flag` BOOL NOT NULL,
    `specified_species` BOOL NOT NULL,
    `hydrogenosome_genetic_code_id` INTEGER NOT NULL,
    `inherited_HGC_flag` BOOL NOT NULL,
    PRIMARY KEY (`tax_id`)
);

CREATE TABLE `division` (
    `division_id` INTEGER NOT NULL,
    `division_cde` CHAR(3) NOT NULL,
    `division_name` VARCHAR NOT NULL,
    `comments` VARCHAR NOT NULL,
    PRIMARY KEY (`division_id`)
);

CREATE TABLE `gencode` (
    `genetic_code_id` INTEGER NOT NULL,
    `abbreviation` VARCHAR,
    `name` VARCHAR NOT NULL,
    `cde` VARCHAR NOT NULL,
    `starts` VARCHAR NOT NULL,
    PRIMARY KEY (`genetic_code_id`)
);

CREATE TABLE `names` (
    `tax_id` INTEGER NOT NULL,
    `name_txt` VARCHAR NOT NULL,
    `unique_name` VARCHAR NOT NULL,
    `name_class` [VARCHAR] NOT NULL,
    `id` SERIAL NOT NULL
);

CREATE TABLE `delnodes` (
    `tax_id` INTEGER NOT NULL,
    PRIMARY KEY (`tax_id`)
);

CREATE TABLE `merged` (
    `old_tax_id` INTEGER NOT NULL,
    `new_tax_id` INTEGER NOT NULL,
    PRIMARY KEY (`new_tax_id`)
);

CREATE TABLE `citations` (
    `cit_id` INTEGER NOT NULL,
    `cit_key` VARCHAR NOT NULL,
    `medline_id` INTEGER NOT NULL,
    `pubmed_id` INTEGER NOT NULL,
    `url` VARCHAR,
    `text` TEXT,
    `taxid_list` VARCHAR,
    PRIMARY KEY (`cit_id`),
    UNIQUE (`medline_id`, `pubmed_id`)
);

CREATE TABLE `typeoftype` (
    `type_name` VARCHAR NOT NULL,
    `synonyms` VARCHAR[] NOT NULL,
    `nomenclature` CHAR(1) NOT NULL,
    `description` TEXT NOT NULL,
    PRIMARY KEY (`type_name`)
);

CREATE TABLE `host` (
    `tax_id` INTEGER NOT NULL,
    `potential_hosts` VARCHAR NOT NULL,
    PRIMARY KEY (`tax_id`)
);

CREATE TABLE `type_material` (
    `tax_id` INTEGER NOT NULL,
    `tax_name` VARCHAR NOT NULL,
    `type` INTEGER NOT NULL,
    `identifier` VARCHAR NOT NULL,
    `id` SERIAL NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE TABLE `rankedlineage` (
    `tax_id` INTEGER NOT NULL,
    `tax_name` VARCHAR NOT NULL,
    `species` VARCHAR NOT NULL,
    `genus` VARCHAR NOT NULL,
    `family` VARCHAR NOT NULL,
    `order` VARCHAR NOT NULL,
    `class` VARCHAR NOT NULL,
    `phylum` VARCHAR NOT NULL,
    `kingdom` VARCHAR NOT NULL,
    `superkingdom` VARCHAR NOT NULL,
    `id` SERIAL NOT NULL
);

CREATE TABLE `fullnamelineage` (
    `tax_id` INTEGER NOT NULL,
    `tax_name` VARCHAR NOT NULL,
    `lineage` VARCHAR NOT NULL,
    PRIMARY KEY (`tax_id`)
);

ALTER TABLE `nodes` ADD FOREIGN KEY (`division_id`) REFERENCES `division`(`division_id`);
ALTER TABLE `nodes` ADD FOREIGN KEY (`genetic_code_id`) REFERENCES `gencode`(`genetic_code_id`);
ALTER TABLE `names` ADD FOREIGN KEY (`tax_id`) REFERENCES `nodes`(`tax_id`);
ALTER TABLE `merged` ADD FOREIGN KEY (`new_tax_id`) REFERENCES `nodes`(`tax_id`);
ALTER TABLE `host` ADD FOREIGN KEY (`tax_id`) REFERENCES `nodes`(`tax_id`);
ALTER TABLE `type_material` ADD FOREIGN KEY (`tax_id`) REFERENCES `nodes`(`tax_id`);
ALTER TABLE `type_material` ADD FOREIGN KEY (`type`) REFERENCES `typeoftype`(`type_name`);
ALTER TABLE `rankedlineage` ADD FOREIGN KEY (`tax_id`) REFERENCES `nodes`(`tax_id`);
ALTER TABLE `fullnamelineage` ADD FOREIGN KEY (`tax_id`) REFERENCES `nodes`(`tax_id`);