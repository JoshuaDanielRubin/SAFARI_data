
./vgan_corrected/dep/vg/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min -q ../euka_dir/euka_db.ry \
                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz  > corrected.gam && echo -e "\n\n"

#./vg_corrected/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min -q ../euka_dir/euka_db.ry \
#                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz > corrected.gam && echo -e "\n\n"

#./vg_uncorrected/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min  \
#                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz > uncorrected.gam && echo -e "\n\n"


#vg stats -a corrected.gam
#echo -e "\n\n"
#vg stats -a uncorrected.gam
