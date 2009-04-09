`produceGEOPlatformFile` <-
function(x.lumi, lib.mapping=NULL, nuIDMode=TRUE, includeAllChipProbe=FALSE, fileName='GEOPlatformFile.txt') {

	if (includeAllChipProbe) {
		chipInfo <- getChipInfo(x.lumi, lib.mapping=lib.mapping)
		chipVersion <- chipInfo$chipVersion[1]
		chipInfo <- getChipInfo(NULL, lib.mapping=lib.mapping, chipVersion=chipVersion, idMapping=TRUE)
	} else {
		chipInfo <- getChipInfo(x.lumi, lib.mapping=lib.mapping, idMapping=TRUE)
		chipVersion <- chipInfo$chipVersion[1]
	}
	idMapping <- chipInfo$idMapping
	chipVersion <- sub("(_V[0-9.]*)_.*$", "\\1", chipVersion)
	organism <- switch(chipInfo$species,
		'Rat'='Rattus norvegicus',
		'Human'="Homo sapiens",
		'Mouse'='Mus musculus')

	platformInfoTitle <- c("Platform_title", "Platform_technology", "Platform_distribution", "Platform_organism", "Platform_manufacturer", "Platform_manufacture_protocol", "Platform_description")
	platformInfo <- c(paste("Illumina expression beadchip", chipVersion), "oligonucleotide beads", "commercial", organism, "Illumina Inc.", "see manufacturer's website", "see manufacturer's website")
	if (nuIDMode) {
		platformInfo[platformInfoTitle == "Platform_title"] <- paste("All Illumina", chipInfo$species,  "expression beadchips based on nuID")
		platformInfo[platformInfoTitle == "Platform_description"] <- paste("This platform pools probe IDs of all Illumina", chipInfo$species,  "expression beadchips based on nuID. nuIDs are defnined based on probe sequence, they are universally unique. Using nuIDs will make data integration and maintenance very easier. See Bioconductor lumi package for more details of nuID annotation.")
	} 
	
	headerTitle <- c("Illumina_ProbeID", "nuID", "Illumina_Gene", "Search_Key", "Accession", "Symbol", "ProbeSequence")
	headerDef <- c("Illumina Probe ID", "nucleotide universal NUcleotide (convertible to and from probe sequence)", "Illumina Gene ID", "Internal id useful for custom design array", "Genbank accession number", "Gene Symbol", "Probe Sequence")
	

	cat('^PLATFORM =', platformInfo[1], '\n', sep='', file=fileName, append=FALSE)
	platformInfoPrint <- paste('!', platformInfoTitle, ' = ', platformInfo, '\n', sep='')
	cat(platformInfoPrint, file=fileName, sep='', append=TRUE)
	mapName <- colnames(idMapping)
	nuID <- idMapping[, 'nuID']
	probeSequence <- id2seq(nuID)
	probe <- idMapping[,grep('probe', mapName, ignore.case=TRUE)]
	gene <- idMapping[,grep('target|gene', mapName, ignore.case=TRUE)]
	accession <- idMapping[,grep('accession', mapName, ignore.case=TRUE)]
	searchKey <- idMapping[,grep('search_key', mapName, ignore.case=TRUE)]
	symbol <- idMapping[,grep('symbol', mapName, ignore.case=TRUE)]
	mappingInfo <- cbind(probe, nuID, gene, searchKey, accession, symbol, probeSequence)

	if (nuIDMode) {
		headerTitle[c(1,2)] <- headerTitle[c(2,1)]		
		headerDef[c(1,2)] <- headerDef[c(2,1)]
		mappingInfo[,c(1,2)] <- mappingInfo[,c(2,1)]
	}
	headerTitle[1] <- 'ID'
	mappingInfo <- rbind(headerTitle, mappingInfo)
	
	platformInfoPrint <- paste('#', headerTitle, ' = ', headerDef, '\n', sep='')
	cat(platformInfoPrint, file=fileName, sep='', append=TRUE)
	cat("!platform_table_begin\n", file=fileName, sep='', append=TRUE)
	write.table(mappingInfo, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)
	cat("!platform_table_end\n", file=fileName, sep='', append=TRUE)
}

