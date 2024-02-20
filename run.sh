#########################
#### Correlation analysis
#########################
output_dir="output/html"
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    echo "Directory $output_dir created."
else
    echo "Directory $output_dir already exists."
fi

file="1-processing"
sequencing="Spatial"
protocole="Spatial" 
docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
        Rscript -e "rmarkdown::render(input = '${file}.Rmd', 
                                      output_file = 'output/html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}'))"

file="1-processing"
sequencing="MinION" 
protocole="scNaUmi_seq" 
docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
        Rscript -e "rmarkdown::render(input = '${file}.Rmd', 
                                      output_file = 'output/html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}'))"

file="1-processing"
sequencing="PromethION"
protocole="scNaUmi_seq" 
docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
        Rscript -e "rmarkdown::render(input = '${file}.Rmd', 
                                      output_file = 'output/html/${sequencing}.html', 
                                      params = list(data_name = '${sequencing}', 
                                      protocole = '${protocole}'))"

file="2-visualization"
docker run --rm -u $(id -u):$(id -g) -v $PWD:/home -w /home genomicpariscentre/sckenver:0.1 \
        Rscript -e "rmarkdown::render(input='${file}.Rmd', output_file='output/html/visualization.html')"


