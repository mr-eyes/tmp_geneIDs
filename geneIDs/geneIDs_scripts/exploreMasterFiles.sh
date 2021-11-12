hgnc_master=$1
ncbi_master=$2
gencode_master=$3

#### Basic exploration of gene catalogs
echo "Basic exploration of gene catalogs"
echo "----------------------------------"
#### HGNC
## Explore
# Explore types of genes
echo "Approved HGNC genes can be calssified as:"
tail -n+2 "$hgnc_master" | awk -F"\t" '{if($6=="Approved")print $4 "|" $5}' | sort | uniq -c
echo "------------------------------------"

#### NCBI
## Explore
## col 10 = type_of_gene, 11 = Symbol_from_nomenclature_authority, 13 = Nomenclature_status
# Explore types of genes:
echo "Current NCBI genes can be calssified as:"
tail -n+2 "$ncbi_master" | awk 'BEGIN{FS="\t";}{print $10}' | sort | uniq -c
echo "NOTE: NCBI IDs labeled as [biological-region] are actually no genes but usually refer to genomic regions with know biological function e.g. trascription regulatory regions"
echo "NOTE: NCBI IDs labeled as [unknown] are usually refering to loci but the exact gene coordinates are not known yet"

# Classification of genes based on naming source (either "O" or "-". The "O" means nomenclature_authority i.e. HGN)
#tail -n+2 "$ncbi_master" | awk 'BEGIN{FS="\t";}{print $13}' | sort | uniq -c ## 19516 - && 42106 O   
#tail -n+2 "$ncbi_master" | awk 'BEGIN{FS="\t";}{if($13=="O")print;}' > Nomenclature_status_O #42106
#grep "HGNC:" Nomenclature_status_O | wc -l     ## 42106
#grep "Ensembl:" Nomenclature_status_O | wc -l  ## 33165
#tail -n+2 "$ncbi_master" | awk 'BEGIN{FS="\t";}{if($13=="O" && $3!=$11)print}' | wc -l # 37
#tail -n+2 "$ncbi_master" | awk 'BEGIN{FS="\t";}{if($13!="O")print;}' > Nomenclature_status_NA #42106
#grep "HGNC:" Nomenclature_status_NA | wc -l     ## 0
#grep "Ensembl:" Nomenclature_status_NA | wc -l  ## 1898
echo "------------------------------------"

#### Gencode_human
## Explore
# Explore types of genes:
echo "Current Gencode genes can be calssified as:"
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $10=="Primary Assembly")print $6}' | sort | uniq -c
echo "------------------------------------"

