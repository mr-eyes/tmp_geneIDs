hgnc_master=$1
ncbi_master=$2
gencode_master=$3
hgnc_his=$4
ncbi_his=$5
gencode_his=$6

#### HGNC
## Simplify the master file for approved HGNC gene IDs to dbXrefs
head -n1 "$hgnc_master" | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$19,$20,$9,$11,$7,$4 "|" $5}' > hgnc_approved.master ## hgnc_id symbol  entrez_id   ensembl_gene_id    alias_symbol    prev_symbol    location   locus_group|locus_type(i.e. gene_type)
cat "$hgnc_master" | awk 'BEGIN{FS=OFS="\t";}{if($19){$19="Entrez:"$19;gsub(",",",Entrez:",$19);}if($6=="Approved")print $1,$2,$19,$20,$9,$11,$7,$4 "|" $5}' >> hgnc_approved.master
hgnc_simple="hgnc_approved.master"

## Stats
echo "Basic HGNC cross reference stats:"
tail -n+2 "$hgnc_simple" | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" approved HGNC ids.", "They are cross referenced with:", "Entrez IDs = "b, "Ensembl IDs = "c}' 
tail -n+2 "$hgnc_simple" | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "HGNC genes with entrez and ensembl ids =",a}'  ## hgnc_id with entrez and ensembl ids = 39578
echo "------------------------------------"


#### NCBI
## Simplify the master file by resolving dbXrefs fields (dbXrefs= $6 in Homo_sapiens.gene_info)
head -n1 "$ncbi_master" | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,"HGNC","Ensembl",$5,$17,$8,$10}' > ncbi_dbXrefs.master ## GeneID  Symbol  HGNC    Ensembl     Synonyms    PrevSymbol   map_location   type_of_gene
tail -n+2 "$ncbi_master" | awk 'BEGIN{FS=OFS="\t"}{if($6!="-")print "Entrez:"$2,$3,$6,$5,$17,$8,$10}' | awk 'BEGIN{FS=OFS="\t"}{ delete vars; split($3,a,"|");for(i in a) { n = index(a[i], ":"); if(n) { x = substr(a[i], n + 1); key = substr(a[i], 1, n - 1); val = substr(a[i], n + 1, length(x)); if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val; } } HGNC = vars["HGNC"]; Ensembl = vars["Ensembl"]; print $1,$2,HGNC,Ensembl,$4,$5,$6,$7; }' >> ncbi_dbXrefs.master
tail -n+2 "$ncbi_master" | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print "Entrez:"$2,$3,"","",$5,$17,$8,$10}' >> ncbi_dbXrefs.master ## $6 = dbXrefs
ncbi_simple="ncbi_dbXrefs.master"

## Stats
echo "Basic NCBI cross reference stats:"
tail -n+2 "$ncbi_simple" | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" current Entrez ids.", "They are cross referenced with:", "HGNC IDs = "b, "and Ensembl IDs = "c}'
tail -n+2 "$ncbi_simple" | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Entrez genes with hgnc and ensembl IDs =",a}'  ## entrez_id with hgnc and ensembl ids= 33165
echo "------------------------------------"


#### Gencode_human
## Simplify the master file for Gencode gene IDs on the Primary Assembly
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,$7,$8,$4,$5,$9,$6}' > gencode_primary.master ## GeneID Symbol  EntrezGene   HGNC    Aliases    PrevSymbols    Location   gene_type
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($7){$7="Entrez:"$7;gsub(",",",Entrez:",$7);}if($3!="Not_in_Gencode" && $10=="Primary Assembly")print $1,$3,$7,$8,$4,$5,$9,$6}' >> gencode_primary.master
genc_simple="gencode_primary.master"

## Stats
echo "Basic GENCODE cross reference stats:"
tail -n+2 "$genc_simple" | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" current Ensembl ids", "They are cross referenced with:", "HGNC IDs = "b, "Entrez IDs = "c}'  
tail -n+2 "$genc_simple" | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Gencode gene with hgnc and entrez  ids =",a}'  ## gencode_id with hgnc and entrez  ids= 24528
echo "------------------------------------"
##################################################################################################################
## Generate 3 pairwise maps from each catalog
tail -n+2 $hgnc_simple | awk 'BEGIN{FS=OFS="\t"}{\
                          if($3){split($3,A,",");for(a in A){print $1,A[a] > "hgnc.hgnc_ncbi.map"}};\
                          if($4){split($4,B,",");for(b in B){print $1,B[b] > "hgnc.hgnc_gencode.map"}};\
                          if($3 && $4){split($3,C,",");split($4,D,",");for(c in C){for(d in D){print C[c],D[d] > "hgnc.ncbi_gencode.map.temp"}}}}'
cat hgnc.ncbi_gencode.map.temp | sort | uniq > hgnc.ncbi_gencode.map
tail -n+2 $ncbi_simple | awk 'BEGIN{FS=OFS="\t"}{\
                          if($3){split($3,A,",");for(a in A){print A[a],$1 > "ncbi.hgnc_ncbi.map"}};\
                          if($4){split($4,B,",");for(b in B){print $1,B[b] > "ncbi.ncbi_gencode.map"}};\
                          if($3 && $4){split($3,C,",");split($4,D,",");for(c in C){for(d in D){print C[c],D[d] > "ncbi.hgnc_gencode.map.temp"}}}}'
cat ncbi.hgnc_gencode.map.temp | sort | uniq > ncbi.hgnc_gencode.map
tail -n+2 $genc_simple | awk 'BEGIN{FS=OFS="\t"}{\
                          if($3){split($3,A,",");for(a in A){print A[a],$1 > "genc.ncbi_gencode.map"}};\
                          if($4){split($4,B,",");for(b in B){print B[b],$1 > "genc.hgnc_gencode.map"}};\
                          if($3 && $4){split($3,C,",");split($4,D,",");for(c in C){for(d in D){print D[d],C[c] > "genc.hgnc_ncbi.map.temp"}}}}'
cat genc.hgnc_ncbi.map.temp | sort | uniq > genc.hgnc_ncbi.map
rm *.map.temp


## count one-to-many relationships
echo "count one-to-many relationships"
echo "-------------------------------"
for map in {hgnc,ncbi,genc}.*.map;do
  sr=$(echo $map | cut -d"." -f1)
  co=$(echo $map | cut -d"." -f2)
  fr=$(echo $map | cut -d"." -f2 | cut -d"_" -f1)
  to=$(echo $map | cut -d"." -f2 | cut -d"_" -f2)
  ##echo $sr $co $fr $to;
  cat $map | awk -F"\t" '{print $1}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print}' > $sr.$co.rep_$fr
  cat $map | awk -F"\t" '{print $2}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print}' > $sr.$co.rep_$to
  #source_simple=$sr"_simple"
  #if [ -s "$sr.$co.rep_$fr" ];then cat $sr.$co.rep_$fr | awk '{print $2}' | grep -Fwf - "${!source_simple}" >  $source_simple.$co.rep.$fr.tab;fi
  #if [ -s "$sr.$co.rep_$to" ];then cat $sr.$co.rep_$to | awk '{print $2}' | grep -Fwf - "${!source_simple}" >  $source_simple.$co.rep.$to.tab;fi
  #cat $sr.$co.rep_$fr | awk '{print $2}' | grep -Fwf - "${!source_simple}" >  $source_simple.$co.rep.$fr.tab.temp
  #cat $sr.$co.rep_$to | awk '{print $2}' | grep -Fwf - "${!source_simple}" >  $source_simple.$co.rep.$to.tab.temp
done
wc -l *.rep_*


## extract one-to-many records 
head -n1 $hgnc_simple > hgnc.many_hgnc;    
cat hgnc.*.rep_hgnc | awk '{print $2}' | grep -Fwf - $hgnc_simple | sort | uniq | sort -t$'\t' -k1,1 >> hgnc.many_hgnc
x=$(cat hgnc.many_hgnc | wc -l); if [ $x -eq 1 ];then rm hgnc.many_hgnc;fi

head -n1 $hgnc_simple > hgnc.many_ncbi;    
cat hgnc.*.rep_ncbi | awk '{print $2}' | grep -Fwf - $hgnc_simple | sort | uniq | sort -t$'\t' -k3,3 >> hgnc.many_ncbi
x=$(cat hgnc.many_ncbi | wc -l); if [ $x -eq 1 ];then rm hgnc.many_ncbi;fi

head -n1 $hgnc_simple > hgnc.many_gencode; 
cat hgnc.*.rep_gencode | awk '{print $2}' | grep -Fwf - $hgnc_simple | sort | uniq | sort -t$'\t' -k4,4 >> hgnc.many_gencode
x=$(cat hgnc.many_gencode | wc -l); if [ $x -eq 1 ];then rm hgnc.many_gencode;fi

head -n1 $ncbi_simple > ncbi.many_hgnc;    
cat ncbi.*.rep_hgnc | awk '{print $2}' | grep -Fwf - $ncbi_simple | sort | uniq | sort -t$'\t' -k3,3  >> ncbi.many_hgnc
x=$(cat ncbi.many_hgnc | wc -l); if [ $x -eq 1 ];then rm ncbi.many_hgnc;fi

head -n1 $ncbi_simple > ncbi.many_ncbi;    
cat ncbi.*.rep_ncbi | awk '{print $2}' | grep -Fwf - $ncbi_simple | sort | uniq | sort -t$'\t' -k1,1 >> ncbi.many_ncbi
x=$(cat ncbi.many_ncbi | wc -l); if [ $x -eq 1 ];then rm ncbi.many_ncbi;fi

head -n1 $ncbi_simple > ncbi.many_gencode; 
cat ncbi.*.rep_gencode | awk '{print $2}' | grep -Fwf - $ncbi_simple | sort | uniq | sort -t$'\t' -k4,4 >> ncbi.many_gencode
x=$(cat ncbi.many_gencode | wc -l); if [ $x -eq 1 ];then rm ncbi.many_gencode;fi

head -n1 $genc_simple > genc.many_hgnc;  
cat genc.*.rep_hgnc | awk '{print $2}' | grep -Fwf - $genc_simple | sort | uniq | sort -t$'\t' -k4,4  >> genc.many_hgnc
x=$(cat genc.many_hgnc | wc -l); if [ $x -eq 1 ];then rm genc.many_hgnc;fi

head -n1 $genc_simple > genc.many_ncbi;
cat genc.*.rep_ncbi | awk '{print $2}' | grep -Fwf - $genc_simple | sort | uniq | sort -t$'\t' -k3,3  >> genc.many_ncbi
x=$(cat genc.many_ncbi | wc -l); if [ $x -eq 1 ];then rm genc.many_ncbi;fi

head -n1 $genc_simple > genc.many_gencode;
cat genc.*.rep_gencode | awk '{print $2}' | grep -Fwf - $genc_simple | sort | uniq | sort -t$'\t' -k1,1  >> genc.many_gencode
x=$(cat genc.many_gencode | wc -l); if [ $x -eq 1 ];then rm genc.many_gencode;fi
echo "-------------------------------"
##################################################################################################################
## cross reference inconsistency 
echo "cross reference inconsistency"
echo "-----------------------------"
## hgnc_ncbi.map
cat {hgnc,ncbi,genc}.hgnc_ncbi.map | sort | uniq -c | sort -k1,1nr | sed 's/ *//' | tr ' ' '\t' > hgnc_ncbi.map_count

cat hgnc_ncbi.map_count | awk -F"\t" '{print $2}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > hgnc_ncbi.dup_hgnc
cat hgnc_ncbi.dup_hgnc | grep -Fwf - *.hgnc_ncbi.map | sed 's/\.hgnc_ncbi\.map:/ /' | tr ' ' '\t' | sort -k2,2 > hgnc_ncbi.dup_hgnc.full
echo hgnc_id HGNC NCBI GENCODE | tr ' ' '\t' > hgnc_ncbi.dup_hgnc.table
cat hgnc_ncbi.dup_hgnc.full | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2][$1])a[$2][$1]=$3;else a[$2][$1]=a[$2][$1]","$3}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> hgnc_ncbi.dup_hgnc.table ## 50
#tail -n+2 hgnc_ncbi.dup_hgnc.table | awk -F"\t" '{if(!$4)print}'
### There are 2 incidents where HGNC assign its id to an NCBI ID while the NCBI assign the HGNC id to another NCBI ID
#tail -n+2 hgnc_ncbi.dup_hgnc.table | awk -F"\t" '{if($4 && $2==$3 && $4!=$2)print}' | wc -l
### There are 48 genes where GENCODE mapped HGNC ID to NCBI ID different from the mapping information in obtained from HGNC or NCBI

cat hgnc_ncbi.map_count | awk -F"\t" '{print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > hgnc_ncbi.dup_ncbi
cat hgnc_ncbi.dup_ncbi | grep -Fwf - *.hgnc_ncbi.map | sed 's/\.hgnc_ncbi\.map:/ /' | tr ' ' '\t' | sort -k3,3 > hgnc_ncbi.dup_ncbi.full
echo ncbi_id HGNC NCBI GENCODE | tr ' ' '\t' > hgnc_ncbi.dup_ncbi.table
cat hgnc_ncbi.dup_ncbi.full | awk 'BEGIN{FS=OFS="\t"}{if(!a["NCBI:"$3][$1])a["NCBI:"$3][$1]=$2;else a["NCBI:"$3][$1]=a["NCBI:"$3][$1]","$2;}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> hgnc_ncbi.dup_ncbi.table ## 23
### 4 of the GENCODE incidents happened because GENCODE assigned an NCBI ID to 2 different HGNC ids

a=$(tail -n+2 hgnc_ncbi.dup_hgnc.table | wc -l)
b=$(cat hgnc_ncbi.dup_hgnc | grep -v -Fwf - hgnc_ncbi.dup_ncbi.table | wc -l)
echo "hgnc_ncbi.map inconsistencies:"$(($a+$b-1))
##########################
## hgnc_gencode.map
cat {hgnc,ncbi,genc}.hgnc_gencode.map | sort | uniq -c | sort -k1,1nr | sed 's/ *//' | tr ' ' '\t' > hgnc_gencode.map_count

cat hgnc_gencode.map_count | awk -F"\t" '{print $2}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > hgnc_gencode.dup_hgnc
cat hgnc_gencode.dup_hgnc | grep -Fwf - *.hgnc_gencode.map | sed 's/\.hgnc_gencode\.map:/ /' | tr ' ' '\t' | sort -k2,2 > hgnc_gencode.dup_hgnc.full
echo hgnc_id HGNC NCBI GENCODE | tr ' ' '\t' > hgnc_gencode.dup_hgnc.table
cat hgnc_gencode.dup_hgnc.full | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2][$1])a[$2][$1]=$3;else a[$2][$1]=a[$2][$1]","$3}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> hgnc_gencode.dup_hgnc.table ## 155

cat hgnc_gencode.map_count | awk -F"\t" '{print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > hgnc_gencode.dup_genc
cat hgnc_gencode.dup_genc | grep -Fwf - *.hgnc_gencode.map | sed 's/\.hgnc_gencode\.map:/ /' | tr ' ' '\t' | sort -k3,3 > hgnc_gencode.dup_genc.full
echo gencode_id HGNC NCBI GENCODE | tr ' ' '\t' > hgnc_gencode.dup_genc.table
cat hgnc_gencode.dup_genc.full | awk 'BEGIN{FS=OFS="\t"}{if(!a[$3][$1])a[$3][$1]=$2;else a[$3][$1]=a[$3][$1]","$2;}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> hgnc_gencode.dup_genc.table ## 63

a=$(tail -n+2 hgnc_gencode.dup_hgnc.table | wc -l)
b=$(cat hgnc_gencode.dup_hgnc | grep -v -Fwf - hgnc_gencode.dup_genc.table | wc -l)
echo "hgnc_gencode.map inconsistencies:"$(($a+$b-1))
##########################
## ncbi_gencode.map
cat {hgnc,ncbi,genc}.ncbi_gencode.map | sort | uniq -c | sort -k1,1nr | sed 's/ *//' | tr ' ' '\t' > ncbi_gencode.map_count

cat ncbi_gencode.map_count | awk -F"\t" '{print $2}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > ncbi_gencode.dup_ncbi
cat ncbi_gencode.dup_ncbi | grep -Fwf - *.ncbi_gencode.map | sed 's/\.ncbi_gencode\.map:/ /' | tr ' ' '\t' | sort -k2,2 > ncbi_gencode.dup_ncbi.full
echo ncbi_id HGNC NCBI GENCODE | tr ' ' '\t' > ncbi_gencode.dup_ncbi.table
cat ncbi_gencode.dup_ncbi.full | awk 'BEGIN{FS=OFS="\t"}{if(!a["NCBI:"$2][$1])a["NCBI:"$2][$1]=$3;else a["NCBI:"$2][$1]=a["NCBI:"$2][$1]","$3;}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> ncbi_gencode.dup_ncbi.table ## 173

cat ncbi_gencode.map_count | awk -F"\t" '{print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > ncbi_gencode.dup_genc
cat ncbi_gencode.dup_genc | grep -Fwf - *.ncbi_gencode.map | sed 's/\.ncbi_gencode\.map:/ /' | tr ' ' '\t' | sort -k3,3 > ncbi_gencode.dup_genc.full
echo gencode_id HGNC NCBI GENCODE | tr ' ' '\t' > ncbi_gencode.dup_genc.table
cat ncbi_gencode.dup_genc.full | awk 'BEGIN{FS=OFS="\t"}{if(!a[$3][$1])a[$3][$1]=$2;else a[$3][$1]=a[$3][$1]","$2;}END{for(i in a){print i,a[i]["hgnc"],a[i]["ncbi"],a[i]["genc"]}}' >> ncbi_gencode.dup_genc.table ## 158


a=$(tail -n+2 ncbi_gencode.dup_ncbi.table | wc -l)
b=$(cat ncbi_gencode.dup_ncbi | grep -v -Fwf - ncbi_gencode.dup_genc.table | wc -l )
echo "ncbi_gencode.map inconsistencies:"$(($a+$b-1))
echo "-----------------------------"
##################################################################################################################
#### 3. Pairwise discrepancies between databases

## NCBI and HGNC
## merge HGNC ids and symbols into the NCBI DB to compare
head -n1 "$ncbi_simple" | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol"}' > NCBI.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 "$hgnc_simple") <(tail -n+2 "$ncbi_simple") >> NCBI.map.hgncExt

## Table 1
echo "## NCBI (Entrez) genes which have HGNC ids and symbols in NCBI different from those in HGNC database:" > hgnc_ncbi.difIDs
echo "Entrez_id,Entrez_Symbol,ncbi-hgnc_id,hgnc-hgnc_id,hgnc-Symbol" >> hgnc_ncbi.difIDs
tail -n+2 NCBI.map.hgncExt | awk 'BEGIN{FS="\t";OFS=","}{if($3!=$9)print $1,$2,($3==""?"-":$3),($9==""?"-":$9),($10==""?"-":$10)}' >> hgnc_ncbi.difIDs

## Table 2
echo "## NCBI (Entrez) genes with the same HGNC ids in NCBI and HGNC database but with different symbols: " > hgnc_ncbi.difSym
echo "Entrez_id,Entrez_Symbol,ncbi-hgnc_id,hgnc-hgnc_id,hgnc-Symbol" >> hgnc_ncbi.difSym
tail -n+2 NCBI.map.hgncExt | awk 'BEGIN{FS="\t";OFS=","}{if($3!="" && $3==$9 && $2!=$10)print $1,$2,$3,$9,($10==""?"-":$10)}' >> hgnc_ncbi.difSym


## Gencode and HGNC
## merge HGNC ids and symbols into Gencode to compare
head -n1 "$genc_simple" | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol","hgnc-hgnc_id2","hgnc-symbol2"}' > gencode.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if(a[$4]=="")a[$4]=$1 FS $2;else a[$4]=a[$4] FS $1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 "$hgnc_simple") <(tail -n+2 "$genc_simple") >> gencode.map.hgncExt

## Table 1
echo "## Gencode genes which have HGNC ids and symbols in Gencode annotation different from those in HGNC database:" > hgnc_gencode.difIDs
echo "Gencode_id,Gencode_Symbol,Gencode-hgnc_id,hgnc-hgnc_id,hgnc-Symbol" >> hgnc_gencode.difIDs
tail -n+2 gencode.map.hgncExt | awk 'BEGIN{FS="\t";OFS=","}{if($4!=$9)print $1,$2,($4==""?"-":$4),($9==""?"-":$9),($10==""?"-":$10)}' >> hgnc_gencode.difIDs

## Table 2
echo "## Gencode genes with the same HGNC ids in Gencode and HGNC database but with different symbols:" > hgnc_gencode.difSym
echo "Gencode_id,Gencode_Symbol,Gencode-hgnc_id,hgnc-hgnc_id,hgnc-Symbol" >> hgnc_gencode.difSym
tail -n+2 gencode.map.hgncExt | awk 'BEGIN{FS="\t";OFS=","}{if($4!="" && $4==$9 && $2!=$10)print $1,$2,$4,$9,($10==""?"-":$10)}' >> hgnc_gencode.difSym

## Gencode and NCBI
## To be implemneted
######################################################
#### cross-referencing to discontinued ID
## HGNC to discontinued NCBI
echo "hgnc_id,discontinued_entrez_id" > hgnc.discont_ncbi
tail -n+2 "$ncbi_his" | awk -F"\t" '{print "Entrez:"$3}' | sort | uniq | grep -Fwf - hgnc.hgnc_ncbi.map | tr '\t' ',' >> hgnc.discont_ncbi
x=$(cat hgnc.discont_ncbi | wc -l); if [ $x -eq 1 ];then rm hgnc.discont_ncbi;fi
## HGNC to discontinued GENCODE
echo "hgnc_id,discontinued_gencode_id" > hgnc.discont_genc
tail -n+2 "$gencode_his" | awk -F"\t" '{print $1}' | sort | uniq | grep -Fwf - hgnc.hgnc_gencode.map | tr '\t' ',' >> hgnc.discont_genc
x=$(cat hgnc.discont_genc | wc -l); if [ $x -eq 1 ];then rm hgnc.discont_genc;fi
## NCBI to discontinued HGNC
echo "discontinued_hgnc_id,entrez_id" > ncbi.discont_hgnc
tail -n+2 "$hgnc_his" | awk -F"\t" '{print $1}' | sort | uniq | grep -Fwf - ncbi.hgnc_ncbi.map | tr '\t' ',' >> ncbi.discont_hgnc
x=$(cat ncbi.discont_hgnc | wc -l); if [ $x -eq 1 ];then rm ncbi.discont_hgnc;fi
## NCBI to discontinued GENCODE
echo "entrez_id,discontinued_gencode_id" > ncbi.discont_genc
tail -n+2 "$gencode_his" | awk -F"\t" '{print $1}' | sort | uniq | grep -Fwf - ncbi.ncbi_gencode.map | tr '\t' ',' >> ncbi.discont_genc
x=$(cat ncbi.discont_genc | wc -l); if [ $x -eq 1 ];then rm ncbi.discont_genc;fi
## GENCODE to discontinued HGNC
echo "discontinued_hgnc_id,gencode_id" > gencode.discont_hgnc
tail -n+2 "$hgnc_his" | awk -F"\t" '{print $1}' | sort | uniq | grep -Fwf - genc.hgnc_gencode.map | tr '\t' ',' >> gencode.discont_hgnc
x=$(cat gencode.discont_hgnc | wc -l); if [ $x -eq 1 ];then rm gencode.discont_hgnc;fi
## GENCODE to discontinued NCBI
echo "discontinued_entrez_id,gencode_id" > gencode.discont_ncbi
tail -n+2 "$ncbi_his" | awk -F"\t" '{print "Entrez:"$3}' | sort | uniq | grep -Fwf - genc.ncbi_gencode.map | tr '\t' ',' >> gencode.discont_ncbi
x=$(cat gencode.discont_ncbi | wc -l); if [ $x -eq 1 ];then rm gencode.discont_ncbi;fi

