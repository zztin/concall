
DIR_INPUT=`realpath $1`
DIR_OUTPUT=`realpath $2` # Directory to store all output in


INS_OR_BB=${DIR_INPUT##*_}
echo 

mkdir -p $DIR_OUTPUT
