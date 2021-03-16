#!/bin/sh
file=$1
groupname=`basename $file .ref`
echo "[ $groupname ] "
grep ^ATOM $file|awk '{if(substr($0, 63,4)=="1.00") {printf(" %d",NR); i++; if(i%10==0) printf("\n")}} END{printf ("\n")}'
