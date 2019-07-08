CID=$(docker create $1)
[ -d output ] && rm -r output
docker cp ${CID}:/benchmark output
docker rm ${CID}

