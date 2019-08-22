-- MySQL Script
-- 
-- Host: localhost    Database: transatlasdb
-- Model: TransAtlasDB		Version: 2.0
-- Function: TransAtlasDB Schema Script
-- 
-- ---------------------------------------------------
-- Server version	5.5.53
/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

-- -----------------------------------------------------
-- Drop all tables if exists (22 tables)
-- -----------------------------------------------------
DROP TABLE IF EXISTS `MapStats`;
DROP TABLE IF EXISTS `GeneStats`;
DROP TABLE IF EXISTS `Metadata`;
DROP TABLE IF EXISTS `VarSummary`;
DROP TABLE IF EXISTS `CommandSyntax`;
DROP TABLE IF EXISTS `VarResult`;
DROP TABLE IF EXISTS `VarAnnotation`;

-- -----------------------------------------------------
-- Table structure for table `MapStats`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `MapStats`;
CREATE TABLE `MapStats` (
	`library_id` int(11) NOT NULL,
	`totalreads` INT(11) NULL DEFAULT NULL,
	`mappedreads` INT(11) NULL DEFAULT NULL,
	`alignmentrate` DOUBLE(5,2) NULL DEFAULT NULL,
	`deletions` INT(11) NULL DEFAULT NULL,
	`insertions` INT(11) NULL DEFAULT NULL,
	`junctions` INT(11) NULL DEFAULT NULL,
	`date` DATE NULL DEFAULT NULL,
	PRIMARY KEY (`library_id`),
	CONSTRAINT `MapStats_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `bird_libraries` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `GeneStats`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `GeneStats`;
CREATE TABLE `GeneStats` (
	`library_id` int(11) NOT NULL,
	`genes` INT(11) NULL DEFAULT NULL,
	`diffexpresstool` VARCHAR(100) NULL DEFAULT NULL,
	`countstool` VARCHAR(100) NULL DEFAULT NULL,
	`date` DATE NULL DEFAULT NULL,
	`countstatus` CHAR(10) NULL,
	`genestatus` CHAR(10) NULL,
	PRIMARY KEY (`library_id`),
	CONSTRAINT `GeneStats_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `MapStats` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `Metadata`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `Metadata`;
CREATE TABLE `Metadata` (
	`library_id` int(11) NOT NULL,
	`refgenome` VARCHAR(100) NULL DEFAULT NULL,
	`annfile` VARCHAR(50) NULL DEFAULT NULL,
	`stranded` VARCHAR(100) NULL DEFAULT NULL,
	`sequencename` TEXT NULL DEFAULT NULL,
	`mappingtool` VARCHAR(100) NULL DEFAULT NULL,
	CONSTRAINT `metadata_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `MapStats` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `VarSummary`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `VarSummary`;
	CREATE TABLE `VarSummary` (
	`library_id` int(11) NOT NULL,
	`totalvariants` INT(11) NULL DEFAULT NULL,
	`totalsnps` INT(11) NULL DEFAULT NULL,
	`totalindels` INT(11) NULL DEFAULT NULL,
	`annversion` VARCHAR(100) NULL DEFAULT NULL,
	`varianttool` VARCHAR(100) NULL DEFAULT NULL,
	`date` DATE NOT NULL, `status` CHAR(10) NULL DEFAULT NULL,
	`nosql` CHAR(10) NULL DEFAULT NULL, PRIMARY KEY (`library_id`),
	CONSTRAINT `varsummary_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `MapStats` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `CommandSyntax`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `CommandSyntax`;
CREATE TABLE `CommandSyntax` (
	`library_id` int(11) NOT NULL,
	`mappingsyntax` TEXT NULL DEFAULT NULL,
	`expressionsyntax` TEXT NULL DEFAULT NULL,
	`countsyntax` TEXT NULL DEFAULT NULL,
	`variantsyntax` TEXT NULL DEFAULT NULL,
	CONSTRAINT `commandsyntax_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `MapStats` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `VarResult`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `VarResult`;
CREATE TABLE `VarResult` (
	`library_id` int(11) NOT NULL,
	`chrom` VARCHAR(100) NOT NULL DEFAULT '',
	`position` INT(11) NOT NULL DEFAULT '0',
	`refallele` VARCHAR(1000) NULL DEFAULT NULL,
	`altallele` VARCHAR(1000) NULL DEFAULT NULL,
	`quality` DOUBLE(20,5) NULL DEFAULT NULL,
	`variantclass` VARCHAR(100) NULL DEFAULT NULL,
	`zygosity` VARCHAR(100) NULL DEFAULT NULL,
	`dbsnpvariant` VARCHAR(100) NULL DEFAULT NULL,
	PRIMARY KEY (`library_id`, `chrom`, `position`),
	CONSTRAINT `varresult_ibfk_1` FOREIGN KEY (`library_id`) REFERENCES `VarSummary` (`library_id`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- Table structure for table `VarAnnotation`
-- -----------------------------------------------------
DROP TABLE IF EXISTS `VarAnnotation`;
CREATE TABLE `VarAnnotation` (
	`library_id` int(11) NOT NULL,
	`chrom` VARCHAR(100) NOT NULL DEFAULT '',
	`position` INT(11) NOT NULL DEFAULT '0',
	`consequence` VARCHAR(100) NOT NULL DEFAULT '',
	`source` VARCHAR(100) NULL DEFAULT NULL,
	`geneid` VARCHAR(1000) NOT NULL DEFAULT '',
	`genename` VARCHAR(1000) NULL DEFAULT NULL,
	`transcript` VARCHAR(250) NULL DEFAULT NULL,
	`feature` VARCHAR(100) NULL DEFAULT NULL,
	`genetype` VARCHAR(250) NULL DEFAULT NULL,
	`proteinposition` VARCHAR(100) NOT NULL DEFAULT '',
	`aachange` VARCHAR(1000) NULL DEFAULT NULL,
	`codonchange` VARCHAR(1000) NULL DEFAULT NULL,
	PRIMARY KEY (`consequence`, `geneid`, `proteinposition`, `library_id`, `chrom`, `position`),
	INDEX `varannotation_indx_genename` (`genename` ASC),
	CONSTRAINT `varannotation_ibfk_1` FOREIGN KEY (`library_id` , `chrom` , `position`) REFERENCES `VarResult` (`library_id` , `chrom` , `position`)
) ENGINE = InnoDB DEFAULT CHARACTER SET = latin1;

-- -----------------------------------------------------
-- View `vw_sampleinfo`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `vw_sampleinfo`;
DROP VIEW IF EXISTS `vw_sampleinfo`;
CREATE VIEW `vw_sampleinfo` AS
	select `a`.`library_id` AS `library_id`, `a`.`species` AS `species`,`a`.`tissue` AS `tissue`,
		`b`.`totalreads` AS `totalreads`, `b`.`mappedreads` AS `mappedreads`, `b`.`alignmentrate` AS `alignmentrate`,
		`c`.`genes` AS `genes`,`d`.`totalvariants` AS `totalvariants`,`d`.`totalsnps` AS `totalsnps`, `d`.`totalindels` AS `totalindels` 
	from (((`bird_libraries` `a` join `MapStats` `b` on((`a`.`library_id` = `b`.`library_id`))) 
		left outer join `GeneStats` `c` on ((`b`.`library_id` = `c`.`library_id`))) 
		left outer join `VarSummary` `d` on ((`a`.`library_id` = `d`.`library_id`)));

-- -----------------------------------------------------
-- View `vw_nosql`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `vw_nosql`;
CREATE VIEW `vw_nosql` AS
	select `a`.`variantclass` AS `variantclass`,`a`.`zygosity` AS `zygosity`,`a`.`dbsnpvariant` AS `dbsnpvariant`,
		`b`.`source` AS `source`,`b`.`consequence` AS `consequence`,`b`.`geneid` AS `geneid`,
		`b`.`genename` AS `genename`,`b`.`transcript` AS `transcript`,`b`.`feature` AS `feature`,
		`b`.`genetype` AS `genetype`,`a`.`refallele` AS `refallele`,`a`.`altallele` AS `altallele`,
		`c`.`tissue` AS `tissue`,`a`.`chrom` AS `chrom`,`b`.`aachange` AS `aachange`,`b`.`codonchange` AS `codonchange`,
		`c`.`species` AS `species`,`a`.`library_id` AS `library_id`,`a`.`quality` AS `quality`,`a`.`position` AS `position`,
		`b`.`proteinposition` AS `proteinposition`
	from ( (`VarResult` `a` join `VarAnnotation` `b` on (( (`a`.`library_id` = `b`.`library_id`) and (`a`.`chrom` = `b`.`chrom`) and (`a`.`position` = `b`.`position`) )) )
		join bird_libraries `c` on ((`a`.`library_id` = `c`.`library_id`)) );

-- -----------------------------------------------------
-- View `vw_vanno`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `vw_vanno`;
CREATE VIEW `vw_vanno` AS
	select `a`.`library_id`, `a`.`chrom` AS `chrom`,`a`.`position` AS `position`,`a`.`refallele` AS `refallele`,
		`a`.`altallele` AS `altallele`,`a`.`variantclass` AS `variantclass`,ifnull(`b`.`consequence`,'-') AS `consequence`,
		`b`.`genename` AS `genename`,`a`.`dbsnpvariant` AS `dbsnpvariant`
	from (`VarResult` `a` left outer join `VarAnnotation` `b` on(((`a`.`library_id` = `b`.`library_id`) and (`a`.`chrom` = `b`.`chrom`) and (`a`.`position` = `b`.`position`))))
		group by `a`.`library_id`, `a`.`chrom`,`a`.`position`,`b`.`consequence`,`b`.`genename`;

-- -----------------------------------------------------
-- View `vw_vvcf`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `vw_vvcf`;
CREATE VIEW `vw_vvcf` AS
	select `a`.`library_id` as `library_id`, `a`.`chrom` AS `chrom`,`a`.`position` AS `position`,
		`a`.`refallele` AS `refallele`,`a`.`altallele` AS `altallele`,`a`.`quality` as `quality`,
		`b`.`consequence` as `consequence`, `b`.`genename` AS `genename`,`b`.`geneid` AS `geneid`,
		`b`.`feature` AS `feature`,`b`.`transcript` AS `transcript`,`b`.`genetype` AS `genetype`,
		`b`.`proteinposition` AS `proteinposition`,`b`.`aachange` AS `aachange`,`b`.`codonchange` AS `codonchange`,
		`a`.`dbsnpvariant` AS `dbsnpvariant`,`a`.`variantclass` AS `variantclass`,`a`.`zygosity` AS `zygosity`,
		`c`.`tissue` AS `tissue`, `c`.`species` AS `species`
	from ((`VarResult` `a` left outer join `VarAnnotation` `b` on (((`a`.`library_id` = `b`.`library_id`) and (`a`.`chrom` = `b`.`chrom`) and (`a`.`position` = `b`.`position`))))
		join `vw_sampleinfo` `c` on ((`a`.`library_id` = `c`.`library_id`))) order by `a`.`library_id`, `a`.`chrom`,`a`.`position`, `b`.`consequence`;

-- -----------------------------------------------------
-- View `vw_seqstats`
-- -----------------------------------------------------
DROP VIEW IF EXISTS `vw_seqstats`;
CREATE VIEW `vw_seqstats` AS
	select `a`.`library_id` AS `library_id`,`a`.`totalreads` AS `totalreads`,`a`.`alignmentrate` AS `alignmentrate`,
		`a`.`genes` AS `genes`,`a`.`totalvariants` AS `totalvariants`,`b`.`mappingtool` AS `mappingtool`,
		`b`.`annfile` AS `annotationfile`,`c`.`date` AS `mapdate`,`d`.`diffexpresstool` AS `diffexpresstool`,
		`d`.`countstool` AS `countstool`,`d`.`date` AS `genedate`,`e`.`varianttool` AS `varianttool`,
		`e`.`annversion` AS `variantannotationtool`,`e`.`date` AS `variantdate`
	from ((((`vw_sampleinfo` `a` join `Metadata` `b` on((`a`.`library_id` = `b`.`library_id`)))
		join `MapStats` `c` on((`a`.`library_id` = `c`.`library_id`)))
		left outer join `GeneStats` `d` on((`a`.`library_id` = `d`.`library_id`)))
		left outer join `VarSummary` `e` on((`a`.`library_id` = `e`.`library_id`))) order by `a`.`library_id`;
-- -----------------------------------------------------
