#### GENCODE genes
echo "Explore the gene symbols ambiguity in GENCODE/Ensembl gene DB"
echo "============================================================="
## Download
wget -O ens_current.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1" >
            
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "external_gene_name" />
        <Attribute name = "external_synonym" />
    </Dataset>
</Query>'
completionStamp=$(tail -n1 ens_current.txt)
if [ "$completionStamp" == "[success]" ]; then
  echo "Ensembl biomart query was done successfuly";
  head -n -1 ens_current.txt > ens_current.txt.temp
  cat ens_current.txt.temp | sort | uniq > ens_current.txt
else echo "Ensembl biomart query is incomplete. Repeat the download step";fi
echo "GeneID" "ensSymbol" "Aliases" | tr ' ' '\t' > ens_current_aggSyn.txt
cat ens_current.txt | awk 'BEGIN{FS=OFS="\t";S="|"}{if(!a[$1FS$2])a[$1FS$2]=$3;else a[$1FS$2]=a[$1FS$2]S$3;}END{for(i in a)print i,a[i]}' >> ens_current_aggSyn.txt

## Find the version of most recent Gencode annotation
curGTF=$(wget -q -O- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/ | \
         grep -o "gencode.v[0-9]\+.chr_patch_hapl_scaff.annotation.gtf.gz" | head -n1)
prefix="gencode.v";
suffix=".chr_patch_hapl_scaff.annotation.gtf.gz";
vcur=${curGTF#$prefix};
vcur=${vcur%$suffix};
#vLast=$(($vcur - 1))
#echo $vcur #$vLast

# Download all assembly report to make assembly map to differentiate 1ry assembly from patches and alternative loci
mkdir -p assemblies_reports && cd assemblies_reports
for sym in M X Y $(seq 1 22);do echo "chr"$sym"|Primary Assembly" | tr '|' '\t';done > assembly_map
wget -q -O- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/ | \
         grep -o "GCA_000001405.[0-9]\+_GRCh3[78].p[0-9]\+" | sort | uniq > assemblies.lst
cat assemblies.lst | while read asm;do
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/${asm}/${asm}_assembly_report.txt
  cat ${asm}_assembly_report.txt | awk 'BEGIN{FS=OFS="\t";}!/#/{print $5,$8}';
  rm ${asm}_assembly_report.txt
done | sort | uniq >> assembly_map
cd ../
  
# Dowanload all previous Gencode GTFs
mkdir -p gencode_gtf && cd gencode_gtf
genFTP="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
declare -A arr; r=2;
for i in 2 2a 3b 3c 3d $(seq 4 ${vcur});do
  arr["$i"]="$r";
  ((r=r+1));
done

extract_geneAnn () {
  i=$1;r=${arr["$i"]};
  gtf=(gencode.v${i}.*.gtf.gz)
  if [ -f "${gtf}" ];then
    #echo $gtf
    #gunzip ${gtf}
    output=${gtf%.gtf.gz}
    zcat ${gtf} | awk -F"\t" '!/#/{if($3=="gene")print $1":"$4"-"$5";"$1";"$9}' | grep -v "_PAR_Y" | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" -v ann_version=${i} -v rank=${r} 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; name = vars["gene_name"]; type = vars["gene_type"]; hgnc = vars["hgnc_id"]; sub(/\..*/,"",id); print id,name,type,hgnc,$1,$2,ann_version,rank; }' | grep -v "^ENSGR" > $output.genes
    echo "GeneID" "Symbol" "gene_type" "HGNC" "Location" "Assembly_type" "gencode_version" "rank" | tr ' ' '\t' > $output.genes.ann
    awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1]=$2;next;}{$6=a[$6];print $0}' ../assemblies_reports/assembly_map $output.genes >> $output.genes.ann
    if [ ! $i -eq ${vcur} ];then rm ${gtf};fi 
    rm $output.genes;
  fi
}

for i in $(seq 20 ${vcur});do #echo $i;
  wget -O gencode.v${i}.ALL.GRCh38.gtf.gz $genFTP/release_${i}/gencode.v${i}.chr_patch_hapl_scaff.annotation.gtf.gz 2>/dev/null; extract_geneAnn "$i";done
for i in $(seq 16 19);do #echo $i;
  wget -O gencode.v${i}.ALL.GRCh37.gtf.gz $genFTP/release_${i}/gencode.v${i}.chr_patch_hapl_scaff.annotation.gtf.gz 2>/dev/null; extract_geneAnn "$i";done
for i in $(seq 5 15);do #echo $i;
  wget -O gencode.v${i}.CHR.GRCh37.gtf.gz $genFTP/release_${i}/gencode.v${i}.annotation.gtf.gz 2>/dev/null; extract_geneAnn "$i";done
wget -O gencode.v4.CHR.GRCh37.gtf.gz $genFTP/release_4/gencode_v4.annotation.GRCh37.gtf.gz 2>/dev/null; extract_geneAnn "4";
wget -O gencode.v3d.CHR.GRCh37.gtf.gz $genFTP/release_3d/gencode.v3d.gtf.gz 2>/dev/null; extract_geneAnn "3d";
wget -O gencode.v3c.CHR.GRCh37.gtf.gz $genFTP/release_3c/gencode.v3c.annotation.GRCh37.gtf.gz 2>/dev/null; extract_geneAnn "3c";
wget -O gencode.v3b.CHR.GRCh37.gtf.gz $genFTP/release_3b/gencode.v3b.annotation.GRCh37.gtf.gz 2>/dev/null; extract_geneAnn "3b";
wget -O gencode.v2a.CHR.GRCh37.gtf.gz $genFTP/release_2/gencode_data.rel2a.gtf.gz 2>/dev/null; extract_geneAnn "2a";
wget -O gencode.v2.CHR.GRCh37.gtf.gz $genFTP/release_2/gencode_data.rel2.gtf.gz 2>/dev/null; extract_geneAnn "2";
wget -O gencode.v1.CHR.GRCh37.gtf.gz $genFTP/release_1/gencode_data.rel1.v2.gtf.gz 2>/dev/null; 

# special processing for release 1
# note: There are 4262 gene IDs in this release without gene symbols
gtf=gencode.v1.CHR.GRCh37.gtf.gz
gunzip ${gtf}
> v1.start; > v1.end;
cat ${gtf%.gz} | awk -F"\t" '!/#/{if($3=="exon")print $1";"$4";"$5";"$9}' | grep -v "_PAR_Y" | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" -v ann_version=1 -v rank=1 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; name = vars["gene_name"]; type = vars["gene_type"]; hgnc = vars["hgnc_id"]; sub(/\..*/,"",id); print id,name,type,hgnc,$1,$2,ann_version,rank >> "v1.start"; print id,$3 >> "v1.end"; }'
cat v1.start | sort -t$'\t' -k1,1 -k6,6  | sort -t$'\t' -k1,1 -u > v1.start.uq
cat v1.end | sort -t$'\t' -k1,1 -k2,2nr  | sort -t$'\t' -k1,1 -u > v1.end.uq
output=${gtf%.gtf.gz}
paste v1.start.uq v1.end.uq | awk 'BEGIN{FS=OFS="\t";}!/#/{print $1,$2,$3,$4,$5":"$6"-"$10,$5,$7,$8;}' | grep -v "^ENSGR" > $output.genes
echo "GeneID" "Symbol" "gene_type" "HGNC" "Location" "Assembly_type" "gencode_version" "rank" | tr ' ' '\t' > $output.genes.ann
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1]=$2;next;}{$6=a[$6];print $0}' ../assemblies_reports/assembly_map $output.genes >> $output.genes.ann
rm gencode.v1.CHR.GRCh37.gtf v1.end v1.start v1.end.uq v1.start.uq $output.genes


cur_Ann=(gencode.v${vcur}.*.genes.ann)
head -n1 $cur_Ann > ../gencode.gene.track
for ann in *.genes.ann;do tail -n+2 $ann;done | sort -t$'\t' -k1,1 -k8,8nr  >> ../gencode.gene.track
for ann in *.genes.ann;do if [ "$ann" != "$cur_Ann" ];then rm $ann;fi;done
cd ../

## update Ensembl gene symbols from current gencode annotation && add all annotation fields from the GTF
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;b[$1]=$3 FS $4 FS $5 FS $6 FS $7;next;}{if(a[$1])print $1,$2,a[$1],$3,b[$1];else print $1,$2,"Not_in_Gencode",$3;}' gencode_gtf/$cur_Ann ens_current_aggSyn.txt > ens_current_aggSyn_genAnn.txt

## generate a uniqe list of previous symbols for each gene ID (check to exclude entries from v1 with no gene symbols)
head -n1 gencode.gene.track > gencode.gene.uniqIDs
tail -n+2 gencode.gene.track | sort -t$'\t' -k1,2 -u | sort -t$'\t' -k1,1 -k8,8nr | awk 'BEGIN{FS=OFS="\t"}{if($2)print}' >> gencode.gene.uniqIDs
cat gencode.gene.uniqIDs |  awk -v vcur=$vcur 'BEGIN{FS=OFS="\t"}{if($7!=vcur)print}' > gencode.gene.uniqIDs.notLast
head -n1 gencode.gene.uniqIDs.notLast > gencode.gene.uniqIDs.notLast.noCur
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1 FS $3]=1;next}{if(!a[$1 FS $2])print}' ens_current_aggSyn_genAnn.txt gencode.gene.uniqIDs.notLast >> gencode.gene.uniqIDs.notLast.noCur
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1 FS $3]=1;next}{if(!a[$1 FS $2])print}' ens_current.txt gencode.gene.uniqIDs.notLast.noCur > gencode.gene.uniqIDs.notLast.noCur.noAlias
echo "GeneID" "PrevSymbols" | tr ' ' '\t' > gencode.gene.aggPrev
tail -n+2 gencode.gene.uniqIDs.notLast.noCur.noAlias | awk 'BEGIN{FS=OFS="\t";S="|"}{if(!a[$1])a[$1]=$2;else a[$1]=a[$1]S$2;}END{for(i in a)print i,a[i]}' >> gencode.gene.aggPrev
awk 'BEGIN{FS=OFS="\t";S="|"}FNR==NR{a[$1]=$2;next}{$4=$4 FS a[$1];print}' gencode.gene.aggPrev ens_current_aggSyn_genAnn.txt > ens_current_aggSyn_aggPrev_genAnn.txt


## Add Entrez IDs
if [ ! -f gencode.v${vcur}.metadata.EntrezGene ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${vcur}/gencode.v${vcur}.metadata.EntrezGene.gz
  gunzip gencode.v${vcur}.metadata.EntrezGene.gz
fi

cur_gtf=(gencode_gtf/gencode.v${vcur}.*.gtf.gz)
zcat $cur_gtf | awk -F"\t" '!/#/{if($3=="transcript")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; trans = vars["transcript_id"]; print trans,id; }' > gencode.trans_to_gene.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{ print a[$1],$2}' gencode.trans_to_gene.map gencode.v${vcur}.metadata.EntrezGene | sort | uniq > gencode.v${vcur}.gen_metadata.EntrezGene
cat gencode.v${vcur}.gen_metadata.EntrezGene | awk 'BEGIN{FS=OFS="\t";S=","}{if(!a[$1])a[$1]=$2;else a[$1]=a[$1]S$2;}END{for(i in a)print i,a[i]}' > gencode.v${vcur}.gen_metadata.EntrezGene.agg
echo "GeneID EntrezGene" | tr ' ' '\t' > gencode.gene_to_entrez.map
cat gencode.v${vcur}.gen_metadata.EntrezGene.agg | awk 'BEGIN{FS="[.\t]";OFS="\t"}{print $1,$3}' | sort | uniq >> gencode.gene_to_entrez.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next;}{$7=a[$1] FS $7;print $0;}' gencode.gene_to_entrez.map ens_current_aggSyn_aggPrev_genAnn.txt > ens_current_aggSyn_aggPrev_genAnn_dbXrefs.txt
gencode_master="ens_current_aggSyn_aggPrev_genAnn_dbXrefs.txt"


## generate a list of symbols for discontinued gene IDs 
cat gencode.gene.uniqIDs.notLast.noCur.noAlias | grep -v ^OTTHUMG | grep -v -Fwf <(tail -n+2 $gencode_master | cut -f1) > gencode.gene.discontinued


#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
echo "GeneID" "Symbol" | tr ' ' '\t' > gencodePrim.ID_to_Current
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $10=="Primary Assembly")print $1,"<"$3">"}' >> gencodePrim.ID_to_Current ## 60665

## Generate a map of current gene IDs to aliases 
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$4}' > gencodePrim.ID_to_EachAlias
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $4!="" && $10=="Primary Assembly")print $1,$4}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> gencodePrim.ID_to_EachAlias

## Generate a map of current gene IDs to previous symbols
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$5}' > gencodePrim.ID_to_EachPrev
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $5!="" && $10=="Primary Assembly")print $1,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> gencodePrim.ID_to_EachPrev

## Generate a map for genes IDs withdrawn to their symbols
echo "GeneID" "Symbol" | tr ' ' '\t' > gencodePrim.ID_to_discontinued
tail -n+2 gencode.gene.discontinued | awk 'BEGIN{FS=OFS="\t";}{if($2!="" && $6=="Primary Assembly")print $1,"<"$2">"}' >> gencodePrim.ID_to_discontinued


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 gencodePrim.ID_to_Current | awk -F"\t" '{print $2}' | sort > gencodePrim.Symbols
## create list of all alias symbols
tail -n+2 gencodePrim.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > gencodePrim.Alias
## create list of all previous symbols
tail -n+2 gencodePrim.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > gencodePrim.Prev
## create list of withdrawn symbols
tail -n+2 gencodePrim.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > gencodePrim.discontinued


#### Basic check
echo "1. Basic check of Ensembl and gencode annotations:"
echo "--------------------------------------------------"
## Problem exploration: Difference between Ensembl and gencode annotations
# 1. Ensembl doesn't have symbols while Gencode does
cat ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3 && $3!="Not_in_Gencode")print}' > ens_current_aggSyn_aggPrev_genAnn_unmatchedSym.txt ##  21606
echo "Number of IDs where Ensembl doesn't have symbols while Gencode does     = " $(tail -n+2 ens_current_aggSyn_aggPrev_genAnn_unmatchedSym.txt | wc -l)       
# 2. Ensembl IDs not found in Gencode annotation:
# Example "ENSG00000274081" belongs to ALT_REF_LOCI. The Ensembl website states that its 1st transcript is part of the Gencode basic annotation
head -n1 ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4}' > ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt
tail -n+2 ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{if($3=="Not_in_Gencode")print $1,$2,$4}' >> ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt ## 123
echo "Number of Ensembl IDs not found in Gencode annotation                   = " $(tail -n+2 ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt | wc -l)       
# 3. Gencode IDs not found in Ensembl annotation
tail -n+2 gencode_gtf/$cur_Ann | awk -F"\t" '{print $1}' | sort > id.gen  ## 67005
tail -n+2 ens_current_aggSyn_aggPrev_genAnn.txt | awk -F"\t" '{print $1}' | sort > id.ens     ## 67128
comm -23 id.gen id.ens > id.gen.sp  ## 0  ##i.e. all gencode IDs present in ensembl
missing_ids=$(cat id.gen.sp | wc -l)
if [ $missing_ids -gt 0 ];then echo "These GENCODE IDS are missing from Ensembl annotation"; cat id.gen.sp;else echo "There are no missing GENCODE IDS from Ensembl annotation";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "Gencode official symbols        = " $(cat gencodePrim.Symbols | wc -l)        ## 60664
echo "Gencode(Ensembl) alias symbols  = " $(cat gencodePrim.Alias | wc -l)          ## 54878
echo "Gencode previous symbols        = " $(cat gencodePrim.Prev | wc -l)           ## 75118
echo "Gencode discontinued symbols    = " $(cat gencodePrim.discontinued | wc -l)   ## 38409
echo "-------------------"


#### check for repeated symbols within the same gene record 
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat gencodePrim.ID_to_Current gencodePrim.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.01.Alias_symbols_matching_current_symbols
cat gencodePrim.ID_to_Current gencodePrim.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.02.Previous_symbols_matching_current_symbols
cat gencodePrim.ID_to_EachAlias gencodePrim.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.03.Previous_symbols_matching_alias_symbols

wc -l Gencodesame.*_matching_*_symbols
echo "Examples for these records (if any):"
cat Gencodesame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced Ensembl IDs"
echo "--------------------------------"
echo "TO BE DONE"
echo "Ongoing effort in temp/ens.trackDiscontinued.sh"
echo "----------------------"

