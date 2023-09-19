
#./vg_corrected/bin/vg minimizer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist -o ../euka_dir/euka_db.min ../euka_dir/euka_db.og -k 13 -w 12 && echo -e "\n\n"
#./vg_corrected/bin/vg rymer -t 60 -p -g ../euka_dir/euka_db.gbwt -d ../euka_dir/euka_db.dist -o ../euka_dir/euka_db.ry ../euka_dir/euka_db.og -k 13 -w 12 && echo -e "\n\n"

#./vgan_corrected/bin/vgan euka -t 60 -fq1 /home/projects2/euka_environments/env_cave_dhigh.fq.gz --euka_dir /home/projects/MAAG/Magpie/euka_dir/ -o euka_corrected_ -a 0.3 -j 0.8

./vgan_uncorrected/bin/vgan euka -t 60 -fq1 /home/projects2/euka_environments/env_cave_damage.pp.fq.gz --euka_dir /home/projects/MAAG/Magpie/euka_dir/ -o euka_uncorrected_


