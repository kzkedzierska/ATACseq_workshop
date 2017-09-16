#!/bash/bin


bedtools intersect -wa -wb -f 0.5 -r -a ./peaks/ENCFF109LQF.bed -b ./peaks/ENCFF146ZCO.bed | sort -k1,1 -k2,2n | bedtools merge > ./peaks/embryoLiver_12.5.bed
bedtools intersect -wa -wb -f 0.5 -r -a ./peaks/ENCFF848NLJ.bed -b ./peaks/ENCFF929LOH.bed | sort -k1,1 -k2,2n | bedtools merge > ./peaks/embryoLiver_0.bed
wc -l ./peaks/*.bed
