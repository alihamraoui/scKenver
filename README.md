# scKeñver
Single-Cell sequencing methodology comparaison workflow

![logo](./sckenver_logo.svg)


scKenver, a term derived from Breton meaning "comparison" or "report". scKenver is designed to aggregate our benchmark analysis of single cell long read methods in a reproducible and automated manner. The purpose of this workflow is to facilitate the replication of results and to ensure maximum transparency in the research process.
# Install
## Docker image

```bash
docker build --rm -t genomicpariscentre/sckenver:0.1 --build-arg R=4.2.2 .
```

