eukadir="/home/projects/MAAG/Magpie/euka_dir/"

#./vgan_corrected/dep/vg/bin/vg giraffe -f  /home/projects2/euka_environments/env_cave_dhigh.fq.gz -t 60 -m $eukadir/euka_db.min -q $eukadir/euka_db.ry \
#                              -b fast -d $eukadir/euka_db.dist -Z $eukadir/euka_db.giraffe.gbz > giraffe_results/corrected.gam && echo -e "\n\n"

./vg_uncorrected/bin/vg giraffe -f  /home/projects2/euka_environments/env_cave_dhigh.fq.gz -t 60 -m $eukadir/euka_db.min \
                              -b fast -d $eukadir/euka_db.dist -Z $eukadir/euka_db.giraffe.gbz  > giraffe_results/uncorrected.gam && echo -e "\n\n"

