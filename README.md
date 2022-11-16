# test neural network FDK
 Quantitatve genetic analysis to evaluate the FDK phenotyping based on a trained neural network
 
Data and R scrips in this repository are used for the analyses described in the manusrctipt entitled 'A Neural Network for Phenotyping Fusarium Damaged Kernels (FDK) and its Impact on Genomic Selection Accuracy'

The genotypic data, in the file 'FHBgeno_2020-21.csv' are arranged with lines in rows and markers in columns. Marker names are in the headers.

The phenotypic data, in the file 'FHBpheno_2020-21.xls' are arranged with column headers in the first line, and each row corresponds to a single experimental unit. The headers are defined as follows:

  GS_trainingset:	A binary variable indicating whether the row of data was included in the genomic selection (GS) model training set.
  GS_validset:	 A binary variable indicating whether the row of data was included in the GS validation set.
  observationUnitDbId:	 The database ID for the experimental unit in the breeding program database.
  studyYear:	 The year when data were collected.
  germplasmName:	 The name of the line being evaluated.
  synonym: Synonyms for the line being evaluated, if they exist.
  name:	 The name of the line being evaluated as specified by the genotyping lab. For lines that were not genotyped, this value is NA.
  studyName:	 The field experiment name.
  plotNumber:	 The plot number within the field experiment
  observationUnitName:	 The name of the plot. This is the combination of experiment name and plot number.
  studyDesign:	 The statistical design of the field experiment.
  replicate:	 The replicate number within the experimental design.
  blockNumber:	 The block number within the experimental design. This is eqivalent to the replicate number in this dataset.
  rowNumber:	 The row position within the experiment where the experimental unit is located.
  colNumber:	 The column position within the experiment where the experimental unit is located.
  DON:	 Observation data for Deoxynivalenol content in units of ppm
  FDK_V:	 Observation data for Fusarium Damaged Kernels (FDK) determined using conventional visual assessments
  FDK_L:	 Observation data for FDK determined using manually labeled images.
  FDK_Lhat:	 Predicted FDK determined using predicted labels based on a neural network.
