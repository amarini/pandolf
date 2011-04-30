#!/bin/bash -f
# usage: ./createLists.sh /pnfs/roma1.infn.it/data/cms/store/user/rahatlou/MC/41xv3DA/

maindir=$1

echo "scanning $maindir"

ls -l $maindir | awk '{print $9}' > datasets.txt
ls -l $maindir | awk '{print "'"$maindir"'" "/" $9}' | xargs -i echo "ls " {} " | grep -v \" 0 \" | awk '{print \"{}/\" \$9}'" > commands.txt 

for i in `ls $maindir`
do 
rm -f files_$i.txt
for j in `ls $maindir/$i`
do
echo dcap://cmsrm-se01.roma1.infn.it$maindir/$i/$j >> files_$i.txt
done
done

#
#N=1 #skip first line
#while read LINE ; do
#    N=$((N+1))
#    echo "Processing $LINE"
#    names[${N}]=$LINE
#done < datasets.txt
#
#
#N=0
#while read LINE ; do
#    N=$((N+1))
#    namescommand[${N}]=$LINE
#    namesnum=${#namescommand}
#done < commands.txt
#
#rm -f finalcommand.sh
#
#for ((i=1;i<$N+1;i++)); do
#    echo ${namescommand[${i}]} ">" files_${names[${i}]}".txt" >> finalcommand.sh
#done
#
#echo "Creating lists. It may take time..."
#
#source finalcommand.sh
#
##rm -f datasets.txt
##rm -f commands.txt
#rm -f finalcommand.sh
#
#echo "LISTS are done "
#
