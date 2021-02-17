#!/bin/bash
sample_id=$1
thresh=$2
nperc="$(seqtk comp result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta | awk '{{Tot+=$2; N+=$9}} END {{print N/Tot*100}}')"
echo "$nperc"
        if (( $(echo "$nperc < $thresh" |bc -l) ));
        then
	echo "Nperc= $nperc, sample {wildcards.sample} will be considered as a Good sample, and labelled accordingly: HighCov"
	sampletag="HighCov"
	else
	echo "Nperc= $nperc, sample {wildcards.sample} will be considered as a Bad sample, and labelled accordingly: LowCov"
	sampletag="LowCov"
fi
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.reads.only.sorted."$sampletag".bam
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.unmapped.reads.only.sorted."$sampletag".bam
mv result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted.bam result/"$sample_id"/Mapping/"$sample_id".alignment.removed.duplicates.mapped.1000trimmed.sorted."$sampletag".bam
sed "s/NC_045512.2/"$sample_id"_"$sampletag"/g" result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta > result/"$sample_id"/variantcall/"$sample_id".consensus."$sampletag".fasta
#rm result/"$sample_id"/variantcall/"$sample_id".old.consensus.fasta
echo "$sample_id = $sampletag with Ns: $nperc % of the genome" > result/"$sample_id"/variantcall/"$sample_id".classification.tab
