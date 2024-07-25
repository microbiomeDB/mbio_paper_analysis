## this bearer token is from my (Danielle's) user account. 
## Replace it w yours, or give me free data
## find yours in devtools when you visit the user dataset space

for file in `ls biom/ | sed 's/\.biom//g'`; do 
	curl -X POST \
		-H "Content-type: multipart/form-data" \
		-H 'Authorization: Bearer eyJhbGciOiJFUzUxMiJ9.eyJzdWIiOiIyMjA2ODU3NjAiLCJpc19ndWVzdCI6ZmFsc2UsImlzcyI6Imh0dHBzOi8vZXVwYXRoZGIub3JnL29hdXRoIiwiYXVkIjoiYXBpQ29tcG9uZW50U2l0ZSIsImF6cCI6ImFwaUNvbXBvbmVudFNpdGUiLCJhdXRoX3RpbWUiOjE3MjA3NTMyMjMsImlhdCI6MTcyMDc1MzIyMywiZXhwIjoxODE1MzYxMjIzLCJwcmVmZXJyZWRfdXNlcm5hbWUiOiJkY2FsbGFuLjIyMDY4NTc2MCIsInNpZ25hdHVyZSI6IjkwMTZiYzljYmY5ZGQ4MjI3NjY5ZTYxOGFlNDRmNTMyIn0.AN0waUzj4izeT6Wo_edFq0-dt48bPscbxZlUq8xlajlhWkr0cK1J_2kaEENSssdmIIzTAc5ZNZBJ_8LFg7I4AiiWAZQzG1b2hTzqgagkrbibGmLzGn9yHr9SlFMBOin4dnybTKaCcTdkjCqsQzHMTwKXCsVpPZxKWiGfiZyQAQZFWKUc' \
		-F meta='{"name":"'$file'-TEST","summary":"'$file': a study from the R package curatedMetagenomicData","description":"'$file': a study from the R package curatedMetagenomicData","datasetType":{"name":"biom","version":"1.0"},"projects":["MicrobiomeDB"],"dependencies":[],"origin":"direct-upload"}' \
		-F file=@biom/$file.biom \
		https://qa.microbiomedb.org/vdi/vdi-datasets; 
done
