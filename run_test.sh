
./vg_corrected/bin/vg minimizer -t 60 -p -g test_graph/test.gbwt -d test_graph/test.dist -o test_graph/test.min test_graph/test.xg -k 7 -w 9 && echo -e "\n\n"

./vg_corrected/bin/vg rymer -t 60 -p -g test_graph/test.gbwt -d test_graph/test.dist -o test_graph/test.ry test_graph/test.xg -k 7 -w 9 && echo -e "\n\n"

./vg_corrected/bin/vg giraffe -f test_graph/test.fq -t 1 -m test_graph/test.min -q test_graph/test.ry \
                              -d test_graph/test.dist -Z test_graph/test.giraffe.gbz > test.gam && echo -e "\n\n"


./vg_uncorrected/bin/vg giraffe -f test_graph/test.fq -t 1 -m test_graph/test.min \
                              -d test_graph/test.dist -Z test_graph/test.giraffe.gbz > test_uncorrected.gam && echo -e "\n\n"

