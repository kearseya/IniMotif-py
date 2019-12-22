# Figures

## Motif Logo



Based on the counts of the kmers scanned from the input DNA sequence file, we could derive a position frequence matrix (PFM) or position weight matrix (PWM). Based on these matrixes, we could generate a logo, where the total height at each position is the information content and the heights of the letters are proportional to their frequencies.



## Hamming Distance figure



This is a scatter plot with x-axis hamming distance, and y-axis kmer frequency/count. Noise, alpha,and colours have been added to the graph to allow for easier distinction between points. The top six most common kmers have been labeled with their reverse compliment in the same colour and there counts next to their respective sequences.



## Position Distribution figure



This is Kernel Density Estimate (KDE) plot showing the relative probabilities of encountering a transcription factor binding sire with the x-axis the relative position along the sequence read.



## Co-occurence Distribution plot 



The relative proportions of encountering a forward and reverse TFBS on the same string.



## Kmer Frequency figure (SELEX only)



Kmer frequency over SELEX-seq rounds is shown in two versions of line graphs with varying y-axis units: f = kmer count/total count and log(f /(1 − f)) with the 6 most common kmers being labeled. A bar graph is also provided to show total number of kmers for each run.
