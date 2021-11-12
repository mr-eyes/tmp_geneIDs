#### Gene ambiguity
## Identify genes with ambiguous alias (the alias is ambiguous if it matches another alias, previous or current gene symbol)
## Gene ambiguity venn diagram
cat entrez.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.01.Current_symbols_matching_other_current_symbols
cat entrez.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat entrez.Alias | uniq) <(cat entrez.Symbols | uniq) > Entrez.03.Alias_symbols_matching_current_symbols
cat entrez.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat entrez.Prev | uniq) <(cat entrez.Symbols | uniq) > Entrez.05.Previous_symbols_matching_current_symbols
comm -12 <(cat entrez.Prev | uniq) <(cat entrez.Alias | uniq)  > Entrez.06.Previous_symbols_matching_alias_symbols
cat entrez.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Symbols | uniq) > Entrez.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Alias | uniq) > Entrez.09.Discontinued_symbols_matching_alias_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Prev | uniq) > Entrez.10.Discontinued_symbols_matching_previous_symbols
cat entrez.replaced | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.11.Replaced_symbols_matching_other_replaced_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Symbols | uniq) > Entrez.12.Replaced_symbols_matching_current_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Alias | uniq) > Entrez.13.Replaced_symbols_matching_alias_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Prev | uniq) > Entrez.14.Replaced_symbols_matching_previous_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.discontinued | uniq) > Entrez.15.Replaced_symbols_matching_discontinued_symbols

echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l Entrez.*_matching_*_symbols

cat entrez.{Symbols,Alias,Prev,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > Entrez.ambiguous_freq.txt
cat Entrez.ambiguous_freq.txt | awk -F"<" '{print "<"$2}' | grep -Fwf - <(cat entrez.ID_to_{Current,EachAlias,EachPrev,discontinued,replaced}) | sort -t$'\t' -k2,2 >  entrez.ambiguous.temp 
head -n1 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"status",$5,$17,"new_symbol_ifReplaced"}' > entrez.complete_withdrawn.temp
tail -n+2 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"Official",$5,$17,"-"}' >> entrez.complete_withdrawn.temp
tail -n+2 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4,$6,$7;}' > withdrawn.temp
awk 'BEGIN{FS=OFS="\t";a["-"]="-";}FNR==NR{a[$2]=$3;next;}{print $1,$2,$3,"-","-",a[$4]}' Homo_sapiens.gene_info withdrawn.temp >> entrez.complete_withdrawn.temp
head -n1 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$2,$3,"status",$5,$17,"new_symbol_ifReplaced"}' > Entrez.ambiguous.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' entrez.complete_withdrawn.temp entrez.ambiguous.temp | sort |  uniq | sort -k1,1 >> Entrez.ambiguous.tab

echo " "
echo "Entrez Gene DB has "$(cat Entrez.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 Entrez.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head Entrez.ambiguous_freq.txt
echo "-------------------------"

