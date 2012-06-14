ls $1/res/*tgz | gawk '{system("tar xzf " $0 )}'
hadd -f mergedToys_$1.root outputToy/hi*root
rm -r outputToy
