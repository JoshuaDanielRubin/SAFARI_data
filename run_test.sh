
./vg_corrected/bin/vg minimizer -t 60 -p -g test_graph/test.gbwt -d test_graph/test.dist -o test_graph/test.min test_graph/test.xg -k 5 -w 5 && echo -e "\n\n"

./vg_corrected/bin/vg rymer -t 60 -p -g test_graph/test.gbwt -d test_graph/test.dist -o test_graph/test.ry test_graph/test.xg -k 5 -w 5 && echo -e "\n\n"

./vg_corrected/bin/vg giraffe -k test_graph/kmers_5.txt -f test_graph/test.fq -t 1 -m test_graph/test.min -q test_graph/test.ry \
                              -d test_graph/test.dist -Z test_graph/test.giraffe.gbz > test.gam && echo -e "\n\n"

vg stats -a test.gam
