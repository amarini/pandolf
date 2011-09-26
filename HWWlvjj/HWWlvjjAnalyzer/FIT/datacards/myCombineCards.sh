#! /bin/bash

DIR=$1

cd $DIR

combineCards.py CMS_hwwlvqq_ele=hwwlvqq_ele.${DIR}.txt CMS_hwwlvqq_mu=hwwlvqq_mu.${DIR}.txt > CMS_hwwlvqq_${DIR}_2channels.txt

cd -
