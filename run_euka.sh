
#./vg_corrected/bin/vg minimizer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist -o ../euka_dir/euka_db.min ../euka_dir/euka_db.og && echo -e "\n\n"

#./vg_corrected/bin/vg rymer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist -o ../euka_dir/euka_db.ry ../euka_dir/euka_db.og && echo -e "\n\n"

#./vgan/bin/vgan euka -t 60 -fq1 test_data/three_dhigh_100.fq.gz --euka_db ../euka_dir/ --outFrag

./vg_corrected/bin/vg giraffe -k ../euka_dir/euka_db.12.kmers -f  test_data/three_dhigh_100.fq.gz -t 60 -m ../euka_dir/euka_db.min -q ../euka_dir/euka_db.ry \
                              -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz > test.gam && echo -e "\n\n"
