#Define function to read the metadata table with the correct formatting
read_metadata = function(file_name){
	read.table(file_name,
						 sep = "\t",
						 header = T,
						 colClasses = c("study" = "character",
														"sample_id" = "character",
														"sample_name" = "character",
														"country" = "character",
														"age" = "numeric",
														"bmi" = "numeric",
														"sex" = "character",
														"t2d" = "character",
														"nafld_case_control" = "character",
														"diagnosis_method" = "character",
														"nafld_score" = "numeric",
														"nafld_alt_score" = "numeric"),
						 row.names = 2,
						 check.names = FALSE,
						 stringsAsFactors = FALSE)
}


read_metaphlan = function(metaphlan_merged_table){
	#This function takes input file path of metaphlan4 and 
	#returns class of matrices on different taxonomic levels
	#kingdom > phylum > class > order > family > genus > species > sgb
	require(dplyr)
	require(textshape)
	raw_table = read.table(metaphlan_merged_table, header = T, sep = "\t", check.names = F)
			## Create taxonomy matrices for different levels
	kingdom = filter(raw_table, grepl("k__", clade_name) & !grepl("p__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	phylum = filter(raw_table, grepl("p__", clade_name) & !grepl("c__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	class = filter(raw_table, grepl("c__", clade_name) & !grepl("o__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	order = filter(raw_table, grepl("o__", clade_name) & !grepl("f__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	family = filter(raw_table, grepl("f__", clade_name) & !grepl("g__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	genus = filter(raw_table, grepl("g__", clade_name) & !grepl("s__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	species = filter(raw_table, grepl("s__", clade_name) & !grepl("t__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
	sgb = filter(raw_table, grepl("t__", clade_name))%>%column_to_rownames("clade_name")%>%as.data.frame()
		
	## Crate taxonomy metadata for different level
		
	kingdom_md = data.frame(type = ifelse(grepl("GB", rownames(kingdom)), "unknown", "known"))%>%cbind(rownames(kingdom))%>%column_to_rownames("rownames(kingdom)")
	phylum_md = data.frame(type = ifelse(grepl("GB", rownames(phylum)), "unknown", "known"))%>%cbind(rownames(phylum))%>%column_to_rownames("rownames(phylum)")
	class_md = data.frame(type = ifelse(grepl("CFGB", rownames(class)), "unknown", "known"))%>%cbind(rownames(class))%>%column_to_rownames("rownames(class)")
	order_md = data.frame(type = ifelse(grepl("OFGB", rownames(order)), "unknown", "known"))%>%cbind(rownames(order))%>%column_to_rownames("rownames(order)")
	family_md = data.frame(type = ifelse(grepl("FGB", rownames(family)), "unknown", "known"))%>%cbind(rownames(family))%>%column_to_rownames("rownames(family)")
	genus_md = data.frame(type = ifelse(grepl("GGB", rownames(genus)), "unknown", "known"))%>%cbind(rownames(genus))%>%column_to_rownames("rownames(genus)")
	species_md = data.frame(type = ifelse(grepl("SGB", rownames(species)), "unknown", "known"))%>%cbind(rownames(species))%>%column_to_rownames("rownames(species)")
	sgb_md = data.frame(type = ifelse(grepl("_SGB.*t__", rownames(sgb)), "unknown", "known"))%>%cbind(rownames(sgb))%>%column_to_rownames("rownames(sgb)")
	
	result = list(
		taxonomy = list(
			kingdom = kingdom,
			phylum = phylum, 
			class = class,
			order = order,
			family = family,
			genus = genus,
			species = species,
			sgb = sgb),
		taxonomy_metadata = list(
			kingdom_metadata = kingdom_md,
			phylum_metadata = phylum_md, 
			class_metadata = class_md,
			order_metadata = order_md,
			family_metadata = family_md,
			genus_metadata = genus_md,
			species_metadata = species_md,
			sgb_metadata = sgb_md)
		)
	
	return(result)
}

metaphlan_pcoa = function(data, metadata , color_column){
	
	jaccard_dist <- vegdist(t(data), method = "jaccard", binary = T)
	jaccard_pcoa <- ecodist::pco(jaccard_dist)
	jaccard_pcoa_df <- data.frame(pcoa1 = jaccard_pcoa$vectors[,1], pcoa2 = jaccard_pcoa$vectors[,2])
	
	jaccard_pcoa_metadata = cbind(jaccard_pcoa_df, metadata)
	
	p = ggplot(jaccard_pcoa_metadata, aes(x=pcoa1, y=pcoa2, color = get(color_column)))+
		geom_point(size=3, alpha = 0.5) +
		labs(x = "pc1",
				 y = "pc2",
				 title = "jaccard pcoa") +
		theme(title = element_text(size = 10))
	return(p)
}

