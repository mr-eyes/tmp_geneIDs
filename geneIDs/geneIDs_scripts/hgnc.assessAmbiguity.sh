#### Gene ambiguity
## Identify genes with ambiguous alias or previous symbol (the alias or previous symbol is ambiguous if it matches another alias, previous or current gene symbol)
## Gene ambiguity venn diagram
cat hgnc.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.01.Current_symbols_matching_other_current_symbols
cat hgnc.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat hgnc.Alias | uniq) <(cat hgnc.Symbols | uniq) > HGNC.03.Alias_symbols_matching_current_symbols
cat hgnc.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat hgnc.Prev | uniq) <(cat hgnc.Symbols | uniq) > HGNC.05.Previous_symbols_matching_current_symbols
comm -12 <(cat hgnc.Prev | uniq) <(cat hgnc.Alias | uniq)  > HGNC.06.Previous_symbols_matching_alias_symbols
cat hgnc.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Symbols | uniq) > HGNC.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Alias | uniq) > HGNC.09.Discontinued_symbols_matching_alias_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Prev | uniq) > HGNC.10.Discontinued_symbols_matching_previous_symbols
cat hgnc.replaced | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.11.Replaced_symbols_matching_other_replaced_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Symbols | uniq) > HGNC.12.Replaced_symbols_matching_current_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Alias | uniq) > HGNC.13.Replaced_symbols_matching_alias_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Prev | uniq) > HGNC.14.Replaced_symbols_matching_previous_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.discontinued | uniq) > HGNC.15.Replaced_symbols_matching_discontinued_symbols

echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l HGNC.*_matching_*_symbols 

cat hgnc.{Symbols,Alias,Prev,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > HGNC.ambiguous_freq.txt
cat HGNC.ambiguous_freq.txt | awk -F"<" '{print "<"$2}' | grep -Fwf - <(cat hgnc.ID_to_{Current,EachAlias,EachPrev,discontinued,replaced}) | sort -t$'\t' -k2,2 >  hgnc.ambiguous.temp 
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11,"new_symbol_ifReplaced"}' > hgnc.complete_withdrawn.temp
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11,"-"}' >> hgnc.complete_withdrawn.temp
tail -n+2 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,$2,"-","-",$4}' >> hgnc.complete_withdrawn.temp
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$1,$2,$6,$9,$11,"new_symbol_ifReplaced"}' > HGNC.ambiguous.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' hgnc.complete_withdrawn.temp hgnc.ambiguous.temp | sort | uniq >> HGNC.ambiguous.tab

echo " "
echo "HGNC has "$(cat HGNC.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 HGNC.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head HGNC.ambiguous_freq.txt
echo "-------------------------"

