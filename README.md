![image](https://user-images.githubusercontent.com/28764332/152460928-d59e1b8c-e3ba-47b0-8557-3e66d8389502.png)

# Liger2LiGer
Nanopore chimera splitting/detection tool

# Dependencies

- build-essentials (CMake v3.10+, make) 
- minimap2 (via PATH or alias)
- samtools v1.9+ (via PATH or alias)

Python3
- matplotlib


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
