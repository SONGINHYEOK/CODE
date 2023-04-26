

import tarfile

tar = tarfile.open("/Users/song-inhyeok/Documents/data/tar_test/tar_total/test_4289558191.tar.gz", "r:gz")

print(tar.getnames())
len(tar.getnames())