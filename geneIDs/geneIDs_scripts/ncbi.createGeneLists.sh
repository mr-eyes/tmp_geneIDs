#### NCBI genes
echo "## Explore the gene symbols ambiguity in Entrez gene DB"
echo "======================================================="
#### Download
if [ ! -f Homo_sapiens.gene_info ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
  gunzip Homo_sapiens.gene_info.gz
fi
if [ ! -f gene_history ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz
  gunzip gene_history.gz
fi

## Update Master files (Part 1):
head -n1 gene_history > human_gene_history 
grep -w ^9606 gene_history | sed 's/~withdrawn//' >> human_gene_history ## The sed command is used to remove "~withdrawn" appended - most probably by mistake - to 17 Discontinued_Symbol 
head -n1 human_gene_history | awk 'BEGIN{FS=OFS="\t";}{print $0,"status","current_GeneID"}' > human_gene_history_track
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$2;next;}!/^#/{if($2=="-"){print $0,"<discontinued>","-";}
                                               else{new_id=$2;while(1){ \
                                                 if(a[new_id]=="-"){print $0,"<discontinued>","-";break;} \
                                                 else if(!a[new_id]){print $0,"<replaced>",new_id;break;} \
                                                 else new_id=a[new_id];}}}' human_gene_history human_gene_history >> human_gene_history_track

#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' > entrez.ID_to_Current
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,"<"$3">"}' >> entrez.ID_to_Current 

## Generate a map of current gene IDs to aliases
head -n1  Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$5}' > entrez.ID_to_EachAlias
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{if($5!="-")print $2,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> entrez.ID_to_EachAlias 

## Generate a map of current gene IDs to previous symbols
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<replaced>")print $7,"<"$4">";}' | sort | uniq > wdEntrez.ID_to_EachPrev
echo "GeneID" "prev_symbol" | tr ' ' '\t' > entrez.ID_to_EachPrev
comm -23 <(sort wdEntrez.ID_to_EachPrev) <(cat entrez.ID_to_Current entrez.ID_to_EachAlias | sort | uniq) >> entrez.ID_to_EachPrev 

## Generate a map for genes IDs withdrawn without approved replacement to their symbols
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4}' > entrez.ID_to_discontinued 
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<discontinued>")print $3,"<"$4">"}' >> entrez.ID_to_discontinued 

## Generate a map for genes IDs withdrawn but replaced by new approved IDs to their symbols
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4}' > entrez.ID_to_replaced
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<replaced>")print $3,"<"$4">"}' >> entrez.ID_to_replaced


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 entrez.ID_to_Current | awk -F"\t" '{print $2}' | sort > entrez.Symbols 
## create list of all alias symbols
tail -n+2 entrez.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > entrez.Alias 
## create list of all previous symbols
tail -n+2 entrez.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > entrez.Prev
## create list of withdrawn symbols without new replacement 
tail -n+2 entrez.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > entrez.discontinued 
## create list of withdrawn symbols with new replacement 
tail -n+2 entrez.ID_to_replaced | awk -F "\t" '{print $2}' | sort > entrez.replaced


## Update Master files (Part 2):
tail -n+2 entrez.ID_to_EachPrev | awk 'BEGIN{FS="[\t<>]";OFS="\t";S="|"}{if(!a[$1])a[$1]=$3;else a[$1]=a[$1]S$3;}END{for(i in a)print i,a[i]}' > entrez.ID_to_EachPrev_aggSyn
awk 'BEGIN{FS=OFS="\t";a["GeneID"]="PrevSymbol";}FNR==NR{a[$1]=$2;next;}{if(a[$2])print $0,a[$2];else print $0,"-";}'  entrez.ID_to_EachPrev_aggSyn Homo_sapiens.gene_info > Homo_sapiens.gene_info.wPrev


echo "1. Basic check of Entrez Ids:"
echo "-----------------------------"
ncbi_IDs=$(tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq | wc -l) ## 61622
ncbi_tot=$(tail -n+2 Homo_sapiens.gene_info | wc -l)
if ((ncbi_tot != ncbi_IDs));then echo "WARNING: The are $ncbi_tot IDs for NCBI genes but only $ncbi_IDs are uniq. The master NCBI dataset has duplicate IDs!!";
else echo "OK! The master NCBI dataset has no duplicate IDs";fi

ncbiHis_IDs=$(tail -n+2 human_gene_history | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq | wc -l)
ncbiHis_tot=$(tail -n+2 human_gene_history | wc -l)
if ((ncbiHis_tot != ncbiHis_IDs));then echo "WARNING: The are $ncbiHis_tot IDs for discontinued NCBI genes but only $ncbiHis_IDs are uniq. There are discontinued NCBI dataset has duplicate IDs!!";
else echo "OK! The withdrawn NBCI dataset has no duplicate IDs";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "Entrez official symbols     = " $(cat entrez.Symbols | wc -l)        ## 63881
echo "Entrez alias symbols        = " $(cat entrez.Alias | wc -l)          ## 71213
echo "Entrez previous symbols     = " $(cat entrez.Prev | wc -l)           ## 15723
echo "Entrez discontinued symbols = " $(cat entrez.discontinued | wc -l)   ## 142208
echo "Entrez replaced symbols     = " $(cat entrez.replaced | wc -l)       ## 21743
echo "-------------------"


## check for repeated symbols in the same gene record
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat entrez.ID_to_Current entrez.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.01.Alias_symbols_matching_current_symbols
cat entrez.ID_to_Current entrez.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.02.Previous_symbols_matching_current_symbols
cat entrez.ID_to_EachAlias entrez.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.03.Previous_symbols_matching_alias_symbols

wc -l Entrezsame.*_matching_*_symbols
echo "Examples for these records (if any):"
cat Entrezsame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced Entrez IDs"
echo "----------------------------------"
echo "No. of withdrawn UIDs replaced by new approved UIDs = "$(cat wdEntrez.ID_to_EachPrev | wc -l)

comm -12 <(sort entrez.ID_to_Current) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev_AsCurrent
echo "No. of withdrawn symbols used as current official symbols = "$(cat wdEntrez.ID_to_EachPrev_AsCurrent | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat wdEntrez.ID_to_EachPrev_AsCurrent | wc -l) "gene IDs"

comm -12 <(sort entrez.ID_to_EachAlias) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as current alias symbols = "$(cat wdEntrez.ID_to_EachPrev_AsAlias | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat wdEntrez.ID_to_EachPrev_AsAlias | wc -l) "gene IDs"

comm -23 <(sort wdEntrez.ID_to_EachPrev) <(cat wdEntrez.ID_to_EachPrev_AsCurrent wdEntrez.ID_to_EachPrev_AsAlias | sort | uniq) > PrevSymbols.missing_From_Current  ## This is the same list we consider "previous symbols"
echo "No. of missing withdrawn symbols that do not show up in the new approved gene records = "$(cat PrevSymbols.missing_From_Current | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat PrevSymbols.missing_From_Current | wc -l) "gene IDs"

grep -v "<LOC.*>" PrevSymbols.missing_From_Current | grep -v "<FLJ.*>" | grep -v "<KIAA.*>" | grep -v "<MGC.*>" | grep -v "<DKFZ[Pp].*>" | grep -vi "<c.*orf.*>" | grep -v "<DJ.*\..*>" | grep -v "<HSA.*>" | grep -v "<PRO.*>" | grep -v "<RP.*-.*\..*>" | grep -v "<AC.*\..*>" | grep -vi "<none>" | grep -vi "<null>" > PrevSymbols.missing_From_Current.knownIDs
echo "No. of missing withdrawn symbols after exclusion of LOC* IDs  = "$(cat PrevSymbols.missing_From_Current.knownIDs | cut -f2 | sort | uniq | wc -l) "symbols for " $(cat PrevSymbols.missing_From_Current.knownIDs | wc -l) "gene IDs"
echo "- Missing symbols that show up as current symbols for other current HGNC genes: "
cat PrevSymbols.missing_From_Current | cut -f2 | grep -Fwf - entrez.ID_to_Current
echo "- Missing symbols that show up as aliases for other current HGNC genes: "
cat PrevSymbols.missing_From_Current | cut -f2 | grep -Fwf - entrez.ID_to_EachAlias
echo "-------------------------"

