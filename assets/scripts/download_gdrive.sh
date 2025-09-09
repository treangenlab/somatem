# run this script in a conda env with conda-forge::gdown installed ; run it from the directory of the script
# downloads the `data/examples/` directory on google drive into the Somatem/examples/ directory

# downloads test data into Somatem/examples/ directory
gdown --folder https://drive.google.com/drive/folders/11ZRpUCRrhdcJarlYdMSEDlCFl3oIz6Bh?usp=sharing -c -O ../../

# note: The `-c` flag to `gdown` will skip already downloaded files. thread: [#99](https://github.com/wkentaro/gdown/issues/99)

# test directory (with 3 text files)
# gdown --folder https://drive.google.com/drive/folders/1vNdGmkeNWSqVpg_Ek6Tdv-kAkqmbuaP_?usp=sharing -O ../../ # test directory
