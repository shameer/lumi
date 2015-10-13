`produceGEOSampleInfoTemplate` <-
function(lumiNormalized, lib.mapping=NULL, fileName='GEOsampleInfo.txt') {	
	
	link <- "http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html"
	templateTitle <- c("Sample_title", "Sample_channel_count","Sample_source_name_ch1","Sample_organism_ch1", "Sample_characteristics_ch1", "Sample_molecule_ch1","Sample_extract_protocol_ch1","Sample_label_ch1", "Sample_label_protocol_ch1", "Sample_hyb_protocol","Sample_scan_protocol","Sample_description","Sample_data_processing","Sample_platform_id", "Sample_supplementary_file")
	if (is(lumiNormalized, 'LumiBatch')) {
		labels <- sampleNames(lumiNormalized)
		chipInfo <- getChipInfo(lumiNormalized, lib.mapping=lib.mapping)
		chipVersion <- chipInfo$chipVersion[1]
		chipVersion <- sub("(_V[0-9.]*)_.*$", "\\1", chipVersion)
		organism <- switch(chipInfo$species,
			'Rat'='Rattus norvegicus',
			'Human'="Homo sapiens",
			'Mouse'='Mus musculus')
		templateContent <- c("","1","",organism,"","total RNA","standard as recommended by illumina","Cy3","standard as recommended by illumina","standard as recommended by illumina","standard as recommended by illumina","","",chipVersion, "none")
	} else if (is(lumiNormalized, 'MethyLumiM')) {
		labels <- sampleNames(lumiNormalized)
		chipVersion <- 'unknown'
		organism <- 'unknown'
		templateContent <- c("","1","",organism,"","genomic DNA","standard as recommended by illumina","Cy3","standard as recommended by illumina","standard as recommended by illumina","standard as recommended by illumina","","",chipVersion,"none")
	} else if (is(lumiNormalized, 'matrix') || is(lumiNormalized, 'ExpressionSet')) {
		if (is(lumiNormalized, 'ExpressionSet'))  lumiNormalized <- exprs(lumiNormalized)
		labels <- colnames(lumiNormalized)
		chipVersion <- 'unknown'
		organism <- 'unknown'
		templateContent <- c("","1","",organism,"","","","Cy3","","standard as recommended by manufacturer","standard as recommended by manufacturer","","", "","none")
	} else {
		cat("The input object should be an object of LumiBatch, MethyLumiM, matrix or other ExpressionSet inherited class!\n")
	}
	
	## add code of parsing processing history 
	preprocessMethod <- ''
	if (is(lumiNormalized, 'LumiBatch')) {
		hh <- getHistory(lumiNormalized)
		comm <- hh$command
	}
	if (is(lumiNormalized, 'MethyLumiM')) {
		hh <- methylumi::getHistory(lumiNormalized)
		comm <- hh$command
	}
	if (is(lumiNormalized, 'LumiBatch')) {
		# grep lumiT method
		lumiT.loc <- grep('lumiT', comm)
		if (length(lumiT.loc) > 0) {
			lumiT.comm <- comm[lumiT.loc]
			method.match <- regexpr('method *= *\"([0-9a-zA-Z]*)\"', lumiT.comm)
			if (method.match == -1) {
				lumiT.method <- 'vst'
			} else {
				end.loc <- method.match + attr(method.match, "match.length") - 2
				method.match <- regexpr('method *= *\"', lumiT.comm)
				start.loc <- method.match + attr(method.match, "match.length")
				lumiT.method <- substring(lumiT.comm, start.loc, end.loc)		
			}
			lumiVersion <- hh$lumiVersion[grep('lumiT', comm)]
			# grep lumiN method
			lumiN.comm <- comm[grep('lumiN', comm)]
			method.match <- regexpr('method *= *\"([0-9a-zA-Z]*)\"', lumiN.comm)
			if (method.match == -1) {
				lumiN.method <- 'quantile'
			} else {
				end.loc <- method.match + attr(method.match, "match.length") - 2
				method.match <- regexpr('method *= *\"', lumiN.comm)
				start.loc <- method.match + attr(method.match, "match.length")
				lumiN.method <- substring(lumiN.comm, start.loc, end.loc)		
			}
			preprocessMethod <- paste('The data was preprocessed by Bioconductor lumi package (version ',lumiVersion, '). It was ', lumiT.method, ' transformed and ', lumiN.method, ' normalized.', sep='')		
		} 		
	} 

	templateContent[templateTitle == "Sample_data_processing"] <- preprocessMethod
	template <- templateTitle
	for (i in seq(labels)) {
		template <- rbind(template, templateContent)
	}
	template <- cbind(c('sampleID', labels), template)
	if (!is.null(fileName)) {
		cat('# For the detailed definition of the column names, please refer to ', link, '\n', sep='', file=fileName)
		write.table(template, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)		
	} else {
		return(template)
	}
}

