# cutout-creator

This repository enables the creation of ERA5 and SARAH weather data Cutouts for PyPSA-Eur: An Open Optimisation Model 
of the European Transmission System. The code is build on top of existing work by https://github.com/fneum.   

## Overview  

The code is built using the [**atlite**](https://atlite.readthedocs.io/en/latest/) Python library and is designed to:  
- Generate separate cutouts for **ERA5** and **SARAH** datasets.  
- Merge and compress the two datasets into a combined cutout.  

## 1. Getting Started    

First, clone the repository:  

```
git clone git@github.com:daniel-rdt/cutout-creator.git
```

## 2. Configure the Settings
Modify the config.yaml file to specify your cutout preferences:

- Cutout: Define the years and parameters for each dataset.
- Merging: Set `merge: true` to merge the ERA5 and SARAH cutouts into a combined and compressed cutout
that will be saved into the folder named `sarah-era5`.
- Existing cutouts: If you already have cutouts, place them into the respective `era5` or `sarah` folders.

## 3. SARAH Cutout Requirements

To generate SARAH cutouts, you must first obtain the SARAH-3 dataset by EUMETSTAT. 
This dataset provides detailed solar radiation data, including:

- Surface Incoming Direct Radiation (SID)
- Surface Incoming Shortwave Radiation (SIS)

It serves as an addition to the ERA5 dataset.

> **_Note:_**
For creating a cutout from this dataset, you need to download large files and your computers memory needs 
> to be able to handle these as well.

### Downloading the SARAH-3 data set

To download the dataset, visit the **EUMETSTAT** website (the link points to the current 3.0 edition):

[Download SARAH-3 Dataset](https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=SARAH_V003)

#### Steps to Download:

1. At the bottom of the page, select the products to include in the cutout:
   - `Surface Incoming Direct Radiation (SID)`
   - `Surface Incoming Shortwave Radiation (SIS)`

   > Select **Instantaneous** time resolution for your chosen year.

2. Add each product to your cart and register on the website.
3. Follow the instructions to activate your account, confirm your order, and wait for the download to be ready.
4. Once notified by email, follow the provided download instructions.
5. Download the ordered files into the `sarah` directory under the following path:
   `sarah/data/sarah-3-<year>`
6. Extract the tar files into the same folder using `tar -xvf *` for Linux systems or use **7zip** for Windows systems.

You are now ready to create cutouts using the SARAH-3 dataset.

## 4. ERA Cutout Requirements
To download and prepare the ERA cutouts, set the configuration option `existing` for ERA to `false`.
Before you can start the processing, you need to complete the two following steps:
- Follow the instructions on the [CDS Copernicus page](https://cds.climate.copernicus.eu/how-to-api) to set up your API access.
This will prompt you to set up an account and place your API key into a local file `$HOME/.cdsapirc`.
- Install the Copernicus Climate Data Store `cdsapi` package in your local environment if not already included.

## 5. Running the code
Finally, to run the code and create the cutouts call:
```
python build_cutout.py
```

> **_Note:_**
Creation of the cutouts as well as merging and compression of the cutouts 
can take up to 6-24 hours and 20-30 GB of memory depending on your local machine.
