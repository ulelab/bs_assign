## bs_assign
Author: aram.amalietti@gmail.com


**Dependencies** (these are the versions the script was developped with, most likely newer version will work but if they don't fall back on these):
```
python=3.7  
pandas=1.2.3  
numpy=1.19.2  
pybedtools=0.8.1  
```
**Usage**:  
  ```
  python3 <path_to_script> <xl_in> <motif1,motif2,...motifn> <prtxn> <fasta> <fai> <window> <len> <cores> <chunk_size> <output_dir> <consensus>  
  ```
  `xl_in` *is a BED file with landmarks around which the motifs are being searched for*  
  `motif1,motif2,...motifn` *is the group of motifs that is searched for around the landmarks, for example AAAA,CCCC,GGGG*  
  `prtxn` *is the path to the file containing relevant positions where the motif need to be in order to be detected, it is an output from another script, PEKA*  
  `fasta` *is the path to the genome in fasta format*  
  `fai` *is the path to the genome index file*  
  `window` *is the distance around the landmarks within which the motifs are being searched for (30 is the usual value)*  
  `len` *is the length of the motifs (4-7 is the usual range)*  
  `cores` *is the number of threads used in the process*  
  `chunk_size` *is the number of landmarks used in one thread (10000 is the usual value)*  
  `output_dir` *is the directory where the results will be saved (make sure it exists)*  
  `consensus` *is a meta motif representing the motif group and will be used in output file names*  
  
**Outputs**:

 BED files representing the binding sites found within the given window around the given landmarks, binding sites are basically found motif merged.  
 there are two merging distances and two corresponding files, one is merged0<file_name> representing merging distance 0 and the other merged30 for 30.
  
  
  
  

