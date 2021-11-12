#### HGNC 
echo "## Explore the gene symbols ambiguity in HGNC"
echo "============================================="
#### Download HGNC dataset (link in the "Statistics & download files" page)
if [ ! -f hgnc_complete_set.txt ];then 
  wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt;fi
if [ ! -f withdrawn.txt ];then
  wget http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt;fi


#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2}' > hgnc.ID_to_Current
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,"<"$2">"}' >> hgnc.ID_to_Current 

## Generate a map of current gene IDs to aliases 
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$9}' > hgnc.ID_to_EachAlias
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($9!="")print $1,$9}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> hgnc.ID_to_EachAlias  

## Generate a map of current gene IDs to previous symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$11}' > hgnc.ID_to_EachPrev
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($11!="")print $1,$11}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> hgnc.ID_to_EachPrev  

## Generate a map for genes IDs withdrawn without approved replacement to their symbols
## Known limitation: We are considering genes replaced by withdrawn genes to be discontinued. There is a possiblity that this new withdrawn gene is also replaced by new approved gene 
##                   but I do not see any example in the current DB (if a gene is replaced, the new genes (column 4) are either "Approved" or "Entry Withdrawn" but not "Merged/Split") 
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' > hgnc.ID_to_discontinued
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Entry Withdrawn")print $1,"<"$3">"}' >> hgnc.ID_to_discontinued;
cat withdrawn.txt | grep -v "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,"<"$3">"}'  >> hgnc.ID_to_discontinued  

## Generate a map for genes IDs withdrawn but replaced by new approved IDs to their symbols
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' > hgnc.ID_to_replaced 
cat withdrawn.txt | grep "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,"<"$3">"}' >> hgnc.ID_to_replaced 


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 hgnc.ID_to_Current | awk -F"\t" '{print $2}' | sort > hgnc.Symbols  
## create list of all alias symbols
tail -n+2 hgnc.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > hgnc.Alias 
## create list of all previous symbols
tail -n+2 hgnc.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > hgnc.Prev
## create list of withdrawn symbols without approved replacement 
tail -n+2 hgnc.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > hgnc.discontinued 
## create list of withdrawn symbols with approved replacement 
tail -n+2 hgnc.ID_to_replaced | awk -F "\t" '{print $2}' | sort > hgnc.replaced


#### Basic check
echo "1. Basic check of HGNC Ids:"
echo "---------------------------"
app_IDs=$(tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
app_tot=$(tail -n+2 hgnc_complete_set.txt | wc -l)
if (( app_tot != app_IDs));then echo "WARNING: The are $app_tot IDs in the approved HGNC but only $app_IDs are uniq. The master HGNC dataset has duplicate IDs!!";
else echo "OK! The master HGNC dataset has no duplicate IDs";fi
nonApp=$(tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{if($6!="Approved")print}' | wc -l)
if (( nonApp > 0));then echo "WARNING: There are $nonApp non-approved records in the master HGNC dataset";
else echo "OK! All records are approved in the master HGNC dataset";fi

wd_IDs=$(tail -n+2 withdrawn.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
wd_tot=$(tail -n+2 withdrawn.txt | wc -l)
if (( wd_tot != wd_IDs));then echo "WARNING: The are $wd_tot IDs in the withdrawn HGNC dataset but only $wd_IDs are uniq. The withdrawn HGNC dataset has duplicate IDs!!";
else echo "OK! The withdrawn HGNC dataset has no duplicate IDs";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "HGNC approved symbols     = " $(cat hgnc.Symbols | wc -l)        ## 42698
echo "HGNC alias symbols        = " $(cat hgnc.Alias | wc -l)          ## 42347
echo "HGNC previous symbols     = " $(cat hgnc.Prev | wc -l)           ## 15180
echo "HGNC discontinued symbols = " $(cat hgnc.discontinued | wc -l)   ## 1826
echo "HGNC replaced symbols     = " $(cat hgnc.replaced | wc -l)       ## 3314
echo "-------------------"


#### check for repeated symbols within the same gene record 
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat hgnc.ID_to_Current hgnc.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.01.Alias_symbols_matching_current_symbols
cat hgnc.ID_to_Current hgnc.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.02.Previous_symbols_matching_current_symbols
cat hgnc.ID_to_EachAlias hgnc.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.03.Previous_symbols_matching_alias_symbols

wc -l HGNCsame.*_matching_*_symbols 
echo "Examples for these records (if any):"
cat HGNCsame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced HGNC IDs"
echo "--------------------------------"
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split"){split($4,a,", ");for(i in a){split(a[i],b,"|");if(b[3]=="Approved")print b[1],"<"$3">";}}}' | sort | uniq > wdHGNC.ID_to_EachPrev ## the terminal sort|uniq is to overcome a bug in the current HGNC report where the HGNC:35188 shows up twice in the same record of the withdrawn HGNC:21128 ID
echo "No. of withdrawn UIDs replaced by new approved UIDs = "$(cat wdHGNC.ID_to_EachPrev | wc -l)

comm -12 <(sort hgnc.ID_to_Current) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsCurrent
echo "No. of withdrawn symbols used as current official symbols = "$(cat wdHGNC.ID_to_EachPrev_AsCurrent | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat wdHGNC.ID_to_EachPrev_AsCurrent | wc -l) "gene IDs"

comm -12 <(sort hgnc.ID_to_EachAlias) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as current alias symbols = "$(cat wdHGNC.ID_to_EachPrev_AsAlias | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat wdHGNC.ID_to_EachPrev_AsAlias | wc -l) "gene IDs"

comm -12 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsPrev
echo "No. of withdrawn symbols reported in previous symbols of current genes = "$(cat wdHGNC.ID_to_EachPrev_AsPrev | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat wdHGNC.ID_to_EachPrev_AsPrev | wc -l) "gene IDs"

comm -23 <(sort wdHGNC.ID_to_EachPrev) <(cat wdHGNC.ID_to_EachPrev_AsCurrent wdHGNC.ID_to_EachPrev_AsAlias wdHGNC.ID_to_EachPrev_AsPrev | sort | uniq) > PrevSymbols.missing_From_Current
echo "No. of missing withdrawn symbols that do not show up in the new approved gene records = "$(cat PrevSymbols.missing_From_Current | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat PrevSymbols.missing_From_Current | wc -l) "gene IDs"
echo "- Missing symbols that show up as current symbols for other current HGNC genes: "
cat PrevSymbols.missing_From_Current | cut -f2 | grep -Fwf - hgnc.ID_to_Current
echo "- Missing symbols that show up as aliases for other current HGNC genes: "
cat PrevSymbols.missing_From_Current | cut -f2 | grep -Fwf - hgnc.ID_to_EachAlias
echo "- Missing symbols that show up as previous symbols for other current HGNC genes: "
cat PrevSymbols.missing_From_Current | cut -f2 | grep -Fwf - hgnc.ID_to_EachPrev

comm -23 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev) > PrevSymbols.missing_From_Withdrawn
echo "No. of previous symbols of current genes that do not show up in the withdrawn symbols = "$(cat PrevSymbols.missing_From_Withdrawn | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat PrevSymbols.missing_From_Withdrawn | wc -l) "gene IDs"
echo "----------------------"

