## find your bearer token in devtools when you visit the user dataset space

for file in `ls biom/ | sed 's/\.biom//g'`; do 
	curl -X POST \
		-H "Content-type: multipart/form-data" \
		-H 'Authorization: [TODO_YOUR_BEARER_TOKEN_HERE]' \
		-F meta='{"name":"'$file'-TEST","summary":"'$file': a study from the R package curatedMetagenomicData","description":"'$file': a study from the R package curatedMetagenomicData","datasetType":{"name":"biom","version":"1.0"},"projects":["MicrobiomeDB"],"dependencies":[],"origin":"direct-upload"}' \
		-F file=@biom/$file.biom \
		https://qa.microbiomedb.org/vdi/vdi-datasets; 
done
