# Comparative analysis of methodologies for detecting extrachromosomal circular DNA
![](figure/figure.png "Overview")
### Description
We created a Python pipeline called ecsim to simulate eccDNA datasets containing both circular and linear DNA from given data. 

### Installation
```
cd $ecsim
pip install .
```

### Usage
```
ecsim --sample SAMPLE_NAME (required) \
      --reference REFERENCE_PATH (required) \
      --thread THREAD_NUMBER (required) \
      --path DATA_PATH (required) \
      --meancov MEAN_DEPTH (default=25) \
      --circular-number CIRCULAR_DNA_NUMBER (default=10000) \
      --linear-number LINEAR_DNA_NUMBER (default=10000) \
      --amp AMPLIFIED_LENGTH (default=5000) \
      --seed RANDOM_SEED (default=None) \
      --simple-ratio SIMPLE_DNA_RATIO (default is calculated from profile) \
      --simple-profile SIMPLE_PROFILE_BEDFILE (containing 5 columns including `chromatin start end length ecID`) \
      --chimeric-profile CHIMERIC_PROFILE_BEDFILE (containing 5 columns including `chromatin start end length ecID`) \
      --ont-model PBSIM2_MODEL_NAME (default='R94') \
      --ont-mean ONT_READ_LENGTH_MEAN (default=3000) \
      --ont-std ONT_READ_LENGTH_STD (default=2500) \
      --sr-platform ILLUMINA_PLATFORM (default='HS25') \
      --sr-mean ILLUMINA_INSERT_LENGTH_MEAN (default=400) \
      --sr-std ILLUMINA_INSERT_LENGTH_STD (default=125) \
      --sr-readlen ILLUMINA_READ_LENGTH (default=150)
```
