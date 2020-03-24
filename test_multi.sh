#5d55c11c10cd48aeba184be1301af2dfcb49ede4 --  test2 always work with this version
./testing.sh test2
echo $?
if [ "$?" == 0 ];
    echo "qsub testing test2 done. successful"
