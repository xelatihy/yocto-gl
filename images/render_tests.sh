mkdir -p tests

for scene in "sh02" "rs01" "ls01" "ps01"
do
    ./bin/yshade --no-ui -r 480 -o tests/$scene.shade.png tests/$scene.obj
done

for scene in "sh03" "rs02" "ls02" "ps02" "cb01"
do
    ./bin/ytrace -r 480 -s 256 --random stratified -o tests/$scene.direct.png tests/$scene.obj
done

for scene in "sh03" "rs02" "ls02" "ps02" "cb01"
do
    ./bin/ytrace -r 480 -s 512 --random stratified -i path -o tests/$scene.path.png tests/$scene.obj
done
