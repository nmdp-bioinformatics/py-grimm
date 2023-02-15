## GRAMM: GRaph bAsed faMilly iMputation.

A tool for imputation haplotypes in families. 

The server is available at: https://gramm.math.biu.ac.il

### Run the code
```
from py-grimm.gram import gram
gram()
```

### Configuration file
The default configuration file is located in 'conf/gram_config.json'.

You can use it, or provide your adjusted configuration file, with the same fields. In this case, provide its absolute path as an argument to gram().


### Input
A CSV file, which includes data of a single family or several families.

A default example file is located in 'source/data/input_file_example.csv'. 
For using your input file, change the value of 'input_file_path' in the configuration file. Use an absolute path.

### Output
1. "results.csv": parent's haplotypes for each family that passed the process successfully.
    Each row is a possible result (max of 100 results for each family).
    If a family has several results, they will appear in descending probability order.
    In addition, in each row, the inheritance (of the children) are described.
2. ".png" files: visualizations of the results: 5 file (or less, if there are less families in the data) for the first 5 families.
    The name of the file is the "FAMCODE" (index) of the family. (e.g. "1.png", "2.png"...).
    If a family has several results, there is a visualization to the first result only.
3. "errors.txt": specifies families, that did not pass the process successfully.
    The file describes the failure reason for each family. 

An output directory is created for every running. Its path is determined by 'output_dir' field in the configuration.
