![image](https://user-images.githubusercontent.com/28764332/152460928-d59e1b8c-e3ba-47b0-8557-3e66d8389502.png)

# Liger2LiGer
Nanopore chimera splitting/detection tool. Automated end-to-end chimera detection and dataset evaluation starting from fastq.

# Dependencies

- build-essentials (CMake v3.10+, make) 
- minimap2 v2.23+ (via PATH or alias)
- samtools v1.9+ (via PATH or alias)

Python3
- matplotlib

# Usage

```
usage: evaluate_chimeras.py [-h] --ref REF --fastq FASTQ [--dry] [--n_threads N_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  --ref REF, -r REF     Path of reference fasta to use for alignment
  --fastq FASTQ, -i FASTQ
                        comma-separated list of fastq file paths to be aligned
  --dry                 echo the commands instead of executing them
  --n_threads N_THREADS
                        How many threads to use for minimap2 alignment
```

Using 64 threads, 6 promethION flowcells (about 1.5TB of sequence data) can be aligned and evaluated in 4.5hr, mainly attributable to recent improvements to minimap2.
![image](https://user-images.githubusercontent.com/28764332/152575429-b908ff7a-f333-4dde-8cc3-caa515c1ac49.png)
(x scale = minutes)

# Output

### Summary

A results_[timestamp].csv with headers:
```
name,n50,n_chimers,n_non_chimers
```
Assuming a unique input filename prefix for each `name` field

### Read IDs

Two files with suffixes:
```
chimeric_reads.txt
non_chimeric_reads.txt
```

With the following format (example):
```
00000131-7ff4-4e70-8b2c-e9f54c3c480a
00000454-bdfa-4a1a-940d-9b6fdf08fb8e
000008c6-0892-471b-a5e1-5072265d40c6
00000bdb-db71-4be4-b21f-326017ceb49d
00000ce6-1638-4fcd-af37-623c8d97c270
00000f22-0a93-4a6f-833b-0d572fb31ea2
```

### Length data

Two files with suffixes:
```
chimer_lengths.txt
non_chimer_lengths.txt
```

With the following format (example):
```
14528
4031
31873
9414
23110
30410
```
where each line contains one read length

### Detailed alignment chains

A list of read IDs and their alignment chains (after splitting).

A chain that has been broken once will appear as:

```
read_id    (A,B),(C,D),
```

where the break exists somewhere between the coordinates B and C in the read sequence. A and D are the start and end of the full alignment.

```
0000275d-8b08-4f25-8c7d-bcc232d4be25	(78,704),(742,1048),(2038,2565),
00008566-66cf-4fe5-8a39-a2cefa93afdf	(88,2080),(2153,2991),(3036,8197),
0001cf78-c05b-4d50-ae53-134fff294ce7	(4691,6563),(32066,36025),
0002e948-b473-48ad-adb0-531ec1add067	(157,2535),(2607,8478),
0003dee8-c3b4-4536-a337-91af6036b473	(41,1514),(1541,7152),
00040b94-98cc-4f1a-90f0-2bcb9bf5bedf	(64,5477),(5519,10855),
```

### Mapping

A `paf` file is stored in the output directory, and it can be split into chimer/non-chimer alignments using the `filter_paf_by_read_name.py` script.


### Plots

Three plots showing the distribution of chimers and non-chimers over a range of read lengths:
- Raw distribution: counts of individual reads
- Coverage distribution: each bin is weighted by its read length (coverage in bp)
- Percent chimeric: percent of reads in each bin that are chimeric


**Example of coverage-weighted distribution:**
![image](https://user-images.githubusercontent.com/28764332/152576219-8e07a297-07bc-44aa-b5bb-5b91d1ad2fc2.png)


# Methods

### Simple case
![image](https://user-images.githubusercontent.com/28764332/152461992-7ccc65a4-bbdc-45e7-9f20-9524eada25e1.png)

### Overlapping case
![image](https://user-images.githubusercontent.com/28764332/152462015-cee614b2-1b7c-4ff3-b76c-3c0cd0a438fe.png)

### Ref contig jump case
![image](https://user-images.githubusercontent.com/28764332/152462086-6d091c88-727d-416f-8b79-14ba9c28f749.png)

### Reversing case
![image](https://user-images.githubusercontent.com/28764332/152462681-20af879e-13f1-4662-bc8c-b2c4e207e545.png)

Option to classify palindromes or split at all reversals is not currently implemented
