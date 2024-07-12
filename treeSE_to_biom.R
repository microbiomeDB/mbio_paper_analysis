writeBIOM <- function(treeSE, filename) {
	
	# currently assumes there is one assay here that we care about
	# and that it is the first one
	# TODO improve that
	assayData <- treeSE@assays@data[[1]]

	rowData <- as.data.frame(SummarizedExperiment::rowData(treeSE))
	# vdi doesnt want the names for some reason :/
	names(rowData) <- NULL 
	# biom constructor does funny stuff if row or col data exist but are empty
	if (ncol(rowData) == 0) rowData <- NULL
	
	colData <- as.data.frame(SummarizedExperiment::colData(treeSE))
	if (ncol(colData) == 0) colData <- NULL

	data_base_type <- ifelse(all(assayData %% 1 == 0), 'int', 'float')

	biom <- biomformat::make_biom(
		assayData, 
		observation_metadata = rowData,
		sample_metadata = colData,
		matrix_element_type = data_base_type
	)

	json <- suppressWarnings(jsonlite::toJSON(biom, always_decimal=TRUE, auto_unbox=TRUE, keep_vec_names=TRUE))
	json <- gsub('"metadata":["', '"metadata":{"taxonomy":["', json, fixed=TRUE)
	json <- gsub(']}', ']}}', json, fixed=TRUE)
	json <- gsub(']]}}', ']]}', json, fixed=TRUE)

	cat(json, file=paste0("biom/", filename,".biom"))
}

setwd('/home/dev')

# getting all of their relative abundance treeSE names
CMD_study_names_list <- curatedMetagenomicData::curatedMetagenomicData(".+relative_abundance")

# keeping most recent
library(dplyr)
studies_matrix <- veupathUtils::strSplit(CMD_study_names_list, '.', 3, c(1,2))
studies_df <- as.data.frame.matrix(studies_matrix)
names(studies_df) <- c('version','study')
studies_df$version <- as.Date(studies_df$version)
keep_studies_df <- studies_df %>% group_by(study) %>% summarize(version=max(version))
keep_studies_names <- paste0(keep_studies_df$version, ".", keep_studies_df$study, ".relative_abundance")

# test w one (for now)
#a_study_name <- keep_studies_names[1]

# getting the actual data now
#a_treeSE <- curatedMetagenomicData::curatedMetagenomicData(a_study_name, dryrun=FALSE)
treeSEs <- lapply(keep_studies_names, curatedMetagenomicData::curatedMetagenomicData, dryrun=FALSE)
treeSEs <- unlist(treeSEs)

# attempting to write them all to biom
#lapply(1:length(a_study_name), function(x){writeBIOM(a_treeSE[[x]], a_study_name[[x]])})
lapply(1:length(keep_studies_names), function(x){writeBIOM(treeSEs[[x]], keep_studies_names[[x]])})
