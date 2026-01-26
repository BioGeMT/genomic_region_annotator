from genomic_region_annotator.download import GTFSpec, download_gtf

if __name__ == "__main__":
    spec = GTFSpec(release=115, species="homo_sapiens", assembly="GRCh38")
    path = download_gtf(spec, "data/annotations")
    print(path)