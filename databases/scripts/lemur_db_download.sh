# lemur : RefSeq v221 bacterial and archaeal genes, and RefSeq v222 fungal genes
curl "https://zenodo.org/records/10802546/files/rv221bacarc-rv222fungi.tar.gz?download=1" --output ../lemur.tar.gz

# unzip the tar file and delete the source
cd ../
tar -xf lemur.tar.gz
rm lemur.tar.gz

# add metadata within the dir
echo "database name: rv221bacarc-rv222fungi.tar.gz" > lemur/db_info
echo "RefSeq v221 bacterial and archaeal genes, and RefSeq v222 fungal genes" >> lemur/db_info
echo "db source: https://zenodo.org/records/10802546/files/rv221bacarc-rv222fungi.tar.gz?download=1" >> lemur/db_info