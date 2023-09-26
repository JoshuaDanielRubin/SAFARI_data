newvg="/home/projects/MAAG/Magpie/Magpie/vg_corrected/bin/vg"

$newvg minimizer -o vgan_corrected/share/vgan/hcfiles/graph.min -t 60 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
                -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og -k 14 -w 12

$newvg rymer -o vgan_corrected/share/vgan/hcfiles/graph.ry -t 60 -d vgan_corrected/share/vgan/hcfiles/graph.dist \
                -g vgan_corrected/share/vgan/hcfiles/graph.gbwt -p vgan_corrected/share/vgan/hcfiles/graph.og -k 14 -w 12
