

#/home/projects/mito_haplotype/vgan/bin/vgan haplocart -t 50 -fq1 test_data/canada_tiny.fq.gz

./vgan_corrected/bin/vgan haplocart -a 0 -t 35 -fq1 test_data/human/B2b3a_damaged.fq --hc-files /home/projects/MAAG/Magpie/Magpie/vgan_corrected/share/vgan/hcfiles -np && echo -e "\n\n"
./vgan_uncorrected/bin/vgan haplocart -t 35 -fq1 test_data/human/B2b3a_damaged.fq --hc-files /home/projects/MAAG/Magpie/Magpie/vgan_corrected/share/vgan/hcfiles -np
