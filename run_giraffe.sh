hcfiles="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/share/vgan/hcfiles"

./vgan_corrected/dep/vg/bin/vg giraffe -f  test_data/human/B2b3a_damaged.fq -t 1 -m $hcfiles/graph.min -q $hcfiles/graph.ry \
                              -b fast -d $hcfiles/graph.dist -Z $hcfiles/graph.giraffe.gbz > corrected.gam && echo -e "\n\n"


#./vg_uncorrected/bin/vg giraffe -f  test_data/human/B2b3a_damaged.fq -t 60 -m $hcfiles/graph.min \
#                              -b fast -d $hcfiles/graph.dist -Z $hcfiles/graph.giraffe.gbz  > uncorrected.gam && echo -e "\n\n"

#./vgan_corrected/dep/vg/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min -q ../euka_dir/euka_db.ry \
#                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz  > corrected.gam && echo -e "\n\n"

#./vg_corrected/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min -q ../euka_dir/euka_db.ry \
#                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz > corrected.gam && echo -e "\n\n"

#./vg_uncorrected/bin/vg giraffe -f  test_data/three_dhigh_100_tinier.fq.gz -t 60 -m ../euka_dir/euka_db.min  \
#                              -b fast -d ../euka_dir/euka_db.dist -Z ../euka_dir/euka_db.giraffe.gbz > uncorrected.gam && echo -e "\n\n"


vg stats -a corrected.gam
echo -e "\n\n"
vg stats -a uncorrected.gam
