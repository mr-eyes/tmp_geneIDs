#### Gene ambiguity
## Identify genes with ambiguous alias (the alias is ambiguous if it matches another alias, previous or current gene symbol)
## Gene ambiguity venn diagram
cat gencodePrim.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.01.Current_symbols_matching_other_current_symbols
cat gencodePrim.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat gencodePrim.Alias | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.03.Alias_symbols_matching_current_symbols
cat gencodePrim.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat gencodePrim.Prev | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.05.Previous_symbols_matching_current_symbols
comm -12 <(cat gencodePrim.Prev | uniq) <(cat gencodePrim.Alias | uniq)  > Gencode.06.Previous_symbols_matching_alias_symbols
cat gencodePrim.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat gencodePrim.discontinued | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat gencodePrim.discontinued | uniq) <(cat gencodePrim.Alias | uniq) > Gencode.09.Discontinued_symbols_matching_alias_symbols
comm -12 <(cat gencodePrim.discontinued | uniq) <(cat gencodePrim.Prev | uniq) > Gencode.10.Discontinued_symbols_matching_previous_symbols


echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l Gencode.*_matching_*_symbols

cat gencodePrim.{Symbols,Alias,Prev,discontinued} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > Gencode.ambiguous_freq.txt
cat Gencode.ambiguous_freq.txt | awk '{print $2}' | grep -Fwf - <(cat gencodePrim.ID_to_{Current,EachAlias,EachPrev,discontinued}) | sort -t$'\t' -k2,2 >  gencodePrim.ambiguous.temp
gencode_master="ens_current_aggSyn_aggPrev_genAnn_dbXrefs.txt"
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,"status",$4,$5,"new_symbol_ifReplaced"}' > gencode.complete_withdrawn.temp
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,"Current",$4,$5,"-"}' > gencode.complete_withdrawn.temp
tail -n+2 gencode.gene.discontinued | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,"discontinued","-","-","-";}' >> gencode.complete_withdrawn.temp
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$1,$3,"status",$4,$5,"new_symbol_ifReplaced"}' > Gencode.ambiguous.tab ## Note that the field "new_symbol_ifReplaced" will be always "-" because we did not map discontinued genes to their new symbols yet.
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' gencode.complete_withdrawn.temp gencodePrim.ambiguous.temp | sort | uniq >> Gencode.ambiguous.tab



echo " "
echo "GENCODE has "$(cat Gencode.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 Gencode.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head Gencode.ambiguous_freq.txt
echo "-------------------------"

