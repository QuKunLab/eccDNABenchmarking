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
      --simple-profile SIMPLE_PROFILE_BEDFILE () \
      --chimeric-profile CHIMERIC_PROFILE_BEDFILE \
      --ont-model PBSIM2_MODEL_NAME \
      --ont-mean ONT_READ_LENGTH_MEAN \
      --ont-std ONT_READ_LENGTH_STD \
      --sr-platform ILLUMINA_PLATFORM \
      --sr-mean ILLUMINA_INSERT_LENGTH_MEAN \
      --sr-std ILLUMINA_INSERT_LENGTH_STD \
      --sr-readlen ILLUMINA_READ_LENGTH
```
