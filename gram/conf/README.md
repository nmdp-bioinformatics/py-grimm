## Configuration file parameters:

| Parameter | Description | Comments |
| --- | --- | --- |
| input_file_path | Path to the input file, which contains the family data. | - |
| alleles_names | Names of the alleles in the data. | Use the default: ["A", "B", "C", "DRB1", "DQB1"] |
| output_dir | Output directory name. | A new directory is created for every running, including the date and time. |
| is_serology_data | Flag - whether the input data is serology. | - |
| races | Dictionary contains families races. The format is: {"CAU", "NAMER"}, for example. If the races are unknown, insert an empty dictionary: {}. | Races only should be used if they are present in grim's configuration file. |
| config_grim_path | Path to grim configuration file. | Use an empty string ("") for using the default grim configuration file (minimal-configuration.json). If you are using an adjusted configuration file, provide its absolute path. |
| build_grim_graph | Flag - whether to run the building of the frequencies graph that grim uses. | Initially, need to be 'true', but after files have been created, can be 'false'. |
