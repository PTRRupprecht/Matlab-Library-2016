This folder contains small scripts that were used to generate a dataset of piano music, downsampled to ca. 7 kHz.

generate_rawData.m contains everything needed to automatically chunk mp3 files in 10-sec pieces (i.e., 63'000 samples per chunk) and save them as mat files.

loadMat_python.py in turn reads some of the mat files and converts them into a dataset.

For more information, please read the related blog entry: https://ptrrupprecht.wordpress.com/2016/09/13/deep-learning-part-iv-deep-dreams-of-music-based-on-dilated-causal-convolutions/
