
import logging
import os
import xarray as xr
import atlite
import yaml
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Sarah
def prepare_cutout_sarah(year, cutout_params):

    cutout = atlite.Cutout(
        path=f"../sarah/europe-{year}-sarah-3.nc",
        module="sarah",
        sarah_dir=f"../sarah/data/sarah-3-{year}",
        chunks={"lat": -1, "time": 100},
        **cutout_params,
    )

    cutout.prepare()

# ERA5
def prepare_cutout_era5(year, cutout_params):

    cutout = atlite.Cutout(
        path=f"../era5/europe-era5-{year}.nc",
        module="era5",
        **cutout_params)

    cutout.prepare(monthly_requests=True)

# Merging
def merge_sarah_with_era5(year):

    era5 = atlite.Cutout(f"../era5/europe-era5-{year}.nc")
    x, y, X, Y = era5.bounds
    sarah = atlite.Cutout(f"../sarah/europe-{year}-sarah-3.nc").sel(x=slice(x,X), y=slice(y,Y))

    cutout = sarah.merge(era5, compat='override')

    os.makedirs("../sarah-era5", exist_ok=True)
    cutout.to_file(f"../sarah-era5/europe-{year}-sarah3-era5.nc")

# Compression
def compress(year):

    ds = xr.open_dataset(f'../sarah-era5/europe-{year}-sarah3-era5.nc')
    comp = dict(zlib=True, complevel=9)  # Maximum compression level
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(f'../sarah-era5/europe-{year}-sarah3-era5-compressed.nc', encoding=encoding)

    #os.remove(f'europe-{year}-sarah3-era5.nc')

if __name__ == "__main__":
    with open('../config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    cutout_params = config["cutout"]
    year = cutout_params["time"]
    cutout_params["x"] = slice(*cutout_params['x'])
    cutout_params["y"] = slice(*cutout_params['y'])

    # Sarah
    if not config["existing_cutouts"]["sarah"]:
        logger.info(f"Preparing SARAH cutout for {year}")
        prepare_cutout_sarah(year, cutout_params)
        logger.info(f"SARAH cutout saved to '/sarah/europe-{year}-sarah-3.nc'")
    else:
        logger.info(f"Cutout for SARAH already prepared and can be imported from '/sarah/europe-{year}-sarah-3.nc'. Please, place input cutout in the required path if not already done.")

    # ERA5
    if not config["existing_cutouts"]["era5"]:
        logger.info(f"Preparing ERA5 cutout for {year}")
        prepare_cutout_era5(year, cutout_params)
        logger.info(f"ERA5 cutout saved to 'era5/europe-era5-{year}.nc'")
    else:
        logger.info(f"Cutout for ERA5 already prepared and can be imported from 'ERA5/europe-era5-{year}.nc'. Please, place input cutout in the required path if not already done.")

    if config["merge"]:
        # Merging
        logger.info(f"Merging SARAH and ERA5 cutout for {year}")
        merge_sarah_with_era5(year)
        logger.info(f"Merged SARAH and ERA5 cutout saved to 'sarah-era5/europe-{year}-sarah3-era5.nc'")

        # Compression
        logger.info(f"Compressing SARAH and ERA5 cutout for {year}")
        compress(year)
        logger.info(f"Compressed SARAH and ERA5 cutout saved to 'sarah-era5/europe-{year}-sarah3-era5-compressed.nc'")