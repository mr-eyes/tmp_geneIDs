#### tracking discontinued GENCODE UIDs
## keep the most recent version of each gene
head -n1 gencode.gene.track | rev | cut -f2- | rev > gencode.gene.annotation
tail -n+2 gencode.gene.track | sort -t$'\t' -k1,1 -u | rev | cut -f2- | rev >> gencode.gene.annotation

## generate history record
cat gencode.gene.annotation | awk -v vcur=$vcur 'BEGIN{FS=OFS="\t"}{if($6!=vcur)print $0}' > gencode.gene.history ## 58776
tail -n+2 gencode.gene.history | awk 'BEGIN{FS=OFS="\t"}{print $6}' | sort | uniq -c | sort -k2,2nr > gencode.gene.history.dist
cat gencode.gene.history | grep -v "^OTTHUMG" > gencode.gene.history.stable  ## 28301
tail -n+2 gencode.gene.history.stable | awk 'BEGIN{FS=OFS="\t"}{print $6}' | sort | uniq -c | sort -k2,2nr > gencode.gene.history.stable.dist
#awk 'BEGIN{OFS="\t"}FNR==NR{a[$2]=$1;next;}{print $2,a[$2],$1}' gencode.gene.history.dist gencode.gene.history.stable.dist > OTTHUMG_dif.txt
mkdir -p his_IDs
tail -n+2 gencode.gene.history.stable | awk 'BEGIN{FS=OFS="\t"}{print $1 > "his_IDs/"$6".ids"}'

## Ensembl Perl API Installation
## https://m.ensembl.org/info/docs/api/api_installation.html
mkdir -p ~/src && cd ~/src
## get bioperl
wget https://github.com/bioperl/bioperl-live/archive/bioperl-release-1-2-3.zip
unzip bioperl-release-1-2-3.zip
mv bioperl-live-bioperl-release-1-2-3 bioperl-1.2.3
rm bioperl-release-1-2-3.zip
## download the Ensembl packages you need
git clone https://github.com/Ensembl/ensembl.git
#git clone https://github.com/Ensembl/ensembl-variation.git
#git clone https://github.com/Ensembl/ensembl-compara.git
#git clone https://github.com/Ensembl/ensembl-funcgen.git
#git clone https://github.com/Ensembl/ensembl-tools.git

echo "export PERL5LIB=$HOME/src/bioperl-1.2.3/:$PERL5LIB" >> ~/.profile
echo "export PERL5LIB=$HOME/src/ensembl/modules:$PERL5LIB" >> ~/.profile
source ~/.profile

wget https://raw.githubusercontent.com/Ensembl/ensembl-tools/release/104/scripts/id_history_converter/IDmapper.pl
chmod +x IDmapper.pl
## Test
#wget https://raw.githubusercontent.com/Ensembl/ensembl-tools/release/104/scripts/id_history_converter/idmapper.in
#./IDmapper.pl -s human -f idmapper.in
#./IDmapper.pl -s human -f gencode.gene.history.ids > gencode.gene.history.raw_idmap
for i in 1 2 2a 3b 3c 3d $(seq 4 $(($vcur - 1)));do
  ./IDmapper.pl -s human -f his_IDs/${i}.ids > his_IDs/${i}.raw_idmap
done
cat his_IDs/*.raw_idmap > gencode.gene.history.raw_idmap


# Assess efficency of ID mapping
cat his_IDs/*.ids | wc -l   ## 28300
cat gencode.gene.history.raw_idmap | grep -v "Old stable ID" | grep -v '^$' | awk -F"." '{print $1}' | uniq | sort > found.ids ## 26643
comm -13 found.ids <(cat his_IDs/*.ids | sort) > missing.ids ## 1657

## Extract the history of Ensembl ID replacement
echo "Old_ID" "New_ID" "ens_Release" "Mapping_score" | tr ' ' '\t' > gencode.gene.history.idmap
cat gencode.gene.history.raw_idmap | sed 's/, /,/g' | sed 's/\./,/' | sed 's/\./,/' | grep -v "^Old" | grep -v "<retired>" | awk -F"," '{if($1!=$3)print $1,$3,$5,$6}' | sort -k1,1 -k4,4dr | sort -k1,2 -u | tr ' ' '\t' >> gencode.gene.history.idmap ## Each old_ID might have one or more new_ID
# annotate
new_tag=$(head -n1 gencode.gene.annotation | awk -v tag="New" 'BEGIN{FS=OFS="\t"}{for(i = 1; i <= NF; ++i)$i=tag"_"$i;sub($1 FS,"");print}')
old_tag=$(head -n1 gencode.gene.annotation | awk -v tag="Old" 'BEGIN{FS=OFS="\t"}{for(i = 1; i <= NF; ++i)$i=tag"_"$i;sub($1 FS,"");print}')
echo "Old_ID" "${old_tag}" "New_ID" "${new_tag}" "ens_Release" "Mapping_score" | tr ' ' '\t' > gencode.gene.history.idmap.ann
notFound=$(echo "-" "-" "-" "-" "-" | tr ' ' '\t')
awk -v nf="${notFound}" 'BEGIN{FS=OFS="\t"}FNR==NR{id=$1;sub($1 FS,"");a[id]=$0;next;}{if(!a[$2])a[$2]=nf;print $1,a[$1],$2,a[$2],$3,$4}' gencode.gene.annotation <(tail -n+2 gencode.gene.history.idmap) >> gencode.gene.history.idmap.ann

## Extract discontinued (or unfound by the tracker) IDs
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next;}{if(!a[$1])print $0;}' gencode.gene.history.idmap gencode.gene.history > gencode.gene.history.unmapped
# add to the annotated IDs
tail -n+2 gencode.gene.history.unmapped | awk 'BEGIN{FS=OFS="\t"}{print $0,"-","-","-","-","-","-","-","-";}' >> gencode.gene.history.idmap.ann

## label discontinued and replaced (Note that the gene will be labled discontinued if it was replaced by gene that was eventually discontinued)
cat gencode.gene.history.idmap.ann | awk 'BEGIN{FS=OFS="\t";}{print $1,$7}' > idmap.temp
cat ens_current_aggSyn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,"current"}' >> idmap.temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next;}{if($7=="New_ID" || $7=="-"){print $0;}
                                               else{new_id=$7;while(1){ \
                                                 if(a[new_id]=="-" || !a[new_id]){print $0,"<discontinued>";break;} \
                                                 else if(a[new_id]=="current"){print $0,"<replaced>";break;} \
                                                 else new_id=a[new_id];}}}' idmap.temp gencode.gene.history.idmap.ann > gencode.gene.history.idmap.details
#grep "replaced"  gencode.gene.history.idmap.details > replaced.temp

