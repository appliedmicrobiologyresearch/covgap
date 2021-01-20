#!/bin/bash
file="$(tail -n +34 $1)"
sampleid=$2
echo "Position  REF	ALT	LABEL	Description	Support" > "$sampleid".secondary_alleles_report.tech
#(head -n 32 result/"$sample_id"/variantcall/"$sample_id".putative.ambdeletions.filteredV3.vcf ; tail -n +33 result/"$sample_id"/variantcall/"$sample_id".putative.ambdeletions.filteredV3.vcf | awk  -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$4,".",$6,$7,$8,$9,$10}') > result/"$sample_id"/variantcall/"$sample_id".putative.ambdeletions.filteredV3corrected.vcf
head -n 31 $1 > "$sampleid".secondary_alleles.vcf
echo "##SECONDARY ALLELE discovery applied on 0 < AF < 0.99" >> "$sampleid".secondary_alleles.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE" >> "$sampleid".secondary_alleles.vcf
while read line
do
CHROM="$(awk '{print $1}' <<< "$line")"
POS="$(awk '{print $2}' <<< "$line")"
ID="$(awk '{print $3}' <<< "$line")"
REF="$(awk '{print $4}' <<< "$line")"
ALT="$(awk '{print $5}' <<< "$line")"
QUAL="$(awk '{print $6}' <<< "$line")"
LABEL="$(awk '{print $7}' <<< "$line")"
INFO="$(awk '{print $8}' <<< "$line")"
QP="$(awk '{print $8}' <<< "$line" | cut -d ';' -f11)"
AF="$(awk '{print $8}' <<< "$line" | cut -d ';' -f2)"
afnumb="$(cut -d '=' -f2 <<< "$AF")"
FORMAT="$(awk '{print $9}' <<< "$line")"
SAMPLE="$(awk '{print $10}' <<< "$line")"
if (( "$(echo "($afnumb < 0.99) && ($afnumb > 0)" |bc -l)" ));
then
qpnumb="$(cut -d '=' -f2 <<< "$QP")"
A="$(cut -d ',' -f1 <<< $qpnumb)"
C="$(cut -d ',' -f2 <<< $qpnumb)"
G="$(cut -d ',' -f3 <<< $qpnumb)"
T="$(cut -d ',' -f4 <<< $qpnumb)"
echo $POS
echo $QP

declare -A QParr=( [A]=$A [C]=$C [G]=$G [T]=$T )

KEYS=$(
for KEY in ${!QParr[@]}; do
  echo "${QParr[$KEY]},$KEY"
done | sort -n -r 
)

tag="dominant"
for i in $KEYS;
do
base="$(cut -d "," -f2 <<< "$i")"
reads="$(cut -d "," -f1 <<< "$i")"
if (( "$(echo "($reads > 0)" |bc -l)" ));
then
echo "$POS	$REF	$base	$LABEL	$tag	$reads" >> "$sampleid".secondary_alleles_report.tech
echo "$CHROM	$POS	$ID	$REF	$base	$QUAL	$LABEL	$INFO	$FORMAT	$SAMPLE" >> "$sampleid".secondary_alleles.vcf
tag="secondary"
else
continue
fi
done

else
continue
fi
done <<< "$file"
