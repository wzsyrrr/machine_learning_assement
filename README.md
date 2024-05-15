# Machine Learning Assessment

This repository contains the code and instructions for the machine learning assessment. Follow the steps below to download the dataset, clean the data, perform analysis, and fit the model.

## Compendium
compendium/ <br/>
├── code    <br/>
│   ├── analysis.R  <br/>
│   └── data_clean.R  <br/>
│   └── model.R  <br/>
├── data_clean  <br/>
│   └── data_clean.rds  <br/>
├── data_raw  <br/>
│   ├── data_raw.txt  <br/>
├── figures  <br/>
│   └── flow_chart.jpeg  <br/>
├── LICENSE  <br/>
├── analysis.sh  <br/>
├── clean.sh  <br/>
├── dwonload.sh  <br/>
├── model.sh  <br/>
└── README.md  <br/>


## Prerequisites

- Access to BlueCrystal 4 (BC4) cluster
- Git installed on your hpc machine

## Steps

### 1. Clone the Repository

First, clone the repository using the following command:

```sh
git clone https://github.com/wzsyrrr/machine_learning_assement.git
```

### 2. Navigate to the Project Directory

Change to the project directory:

```sh
cd machine_learning_assement
```

### 3. Download the Dataset

Run the provided shell script to download the dataset from Google Drive we set up. The data will be saved in the data_raw directory.

```sh
sbatch download.sh
```

### 4. Clean the Data

Use the SLURM job scheduler to clean the data. The cleaned data will be saved in the data_clean directory.

```sh
sbatch clean.sh
```

### 5. Perform Data Analysis

After cleaning the data, run the analysis on the cleaned dataset:

```sh
sbatch analysis.sh
```

### 6. Fit the Model

Once the analysis is complete, fit the model based on the observations from the analysis:

```sh
sbatch model.sh
```

the code will plot results in to the figures file.