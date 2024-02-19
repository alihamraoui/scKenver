#########################
#### Correlation analysis
#########################

sequencing="PromethION" # MinION, PromethION, Spatial
protocole="scNaUmi_seq"  # scNaUmi_seq, Spattial

#sequencing    protocole
#MinION          scNaUmi_seq
#PromethION      scNaUmi_seq
#Spatial         Spatial

file="1-processing"

#docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
#        Rscript -e "rmarkdown::render(input = '${file}.Rmd', 
        #                              output_file = 'output/html/${sequencing}.html', 
        #                              params = list(data_name = '${sequencing}', 
        #                              protocole = '${protocole}'))"

file="2-visualization"

docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
        Rscript -e "rmarkdown::render(input='${file}.Rmd', output_file='output/html/visualization.html')"


