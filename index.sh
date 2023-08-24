newvg="/home/projects/MAAG/Magpie/Magpie/vg_corrected/bin/vg"

$newvg minimizer -o vgan_corrected/share/vgan/hcfiles/graph.min -t 60 -k 12 -w 10 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
                -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og

$newvg rymer -o vgan_corrected/share/vgan/hcfiles/graph.ry -t 60 -k 12 -w 10 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
                -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og
