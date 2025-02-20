
import logging
import os
import xarray as xr
import atlite
import yaml
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Merging
def merge_sarah_with_era5(year):

    era5 = atlite.Cutout(f"../era5/europe-era5-{year}.nc")
    x, y, X, Y = era5.bounds
    sarah = atlite.Cutout(f"../sarah/europe-{year}-sarah-3.nc").sel(x=slice(x,X), y=slice(y,Y))

    cutout = sarah.merge(era5, compat='override')

    os.makedirs("../sarah-era5", exist_ok=True)
    cutout.to_file(snakemake.output.merged)

# Compression
def compress(year):

    ds = xr.open_dataset(f'../sarah-era5/europe-{year}-sarah3-era5.nc')
    comp = dict(zlib=True, complevel=9)  # Maximum compression level
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(snakemake.output.merged_compressed, encoding=encoding)

    #os.remove(f'europe-{year}-sarah3-era5.nc')

if __name__ == "__main__":
    if "snakemake" in globals():
        config = snakemake.config
    else:
        with open('../config.yaml', 'r') as file:
            config = yaml.safe_load(file)

    cutout_params = config["cutout"]
    year = snakemake.wildcards.year if "snakemake" in globals() else cutout_params["time"]
    cutout_params["x"] = slice(*cutout_params['x'])
    cutout_params["y"] = slice(*cutout_params['y'])

    # Merging
    logger.info(f"Merging SARAH and ERA5 cutout for {year}")
    merge_sarah_with_era5(year)
    logger.info(f"Merged SARAH and ERA5 cutout saved to 'sarah-era5/europe-{year}-sarah3-era5.nc'")

    # Compression
    logger.info(f"Compressing SARAH and ERA5 cutout for {year}")
    compress(year)
    logger.info(f"Compressed SARAH and ERA5 cutout saved to 'sarah-era5/europe-{year}-sarah3-era5-compressed.nc'")