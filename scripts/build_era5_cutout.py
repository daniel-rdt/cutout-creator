
import logging
import os
import xarray as xr
import atlite
import yaml
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ERA5
def prepare_cutout_era5(year, cutout_params):

    cutout = atlite.Cutout(
        path=f"snakemake.output",
        module="era5",
        **cutout_params)

    cutout.prepare()

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

    # ERA5
    if not config["existing_cutouts"]["era5"]:
        logger.info(f"Preparing ERA5 cutout for {year}")
        prepare_cutout_era5(year, cutout_params)
        logger.info(f"ERA5 cutout saved to 'era5/europe-era5-{year}.nc'")
    else:
        logger.info(f"Cutout for ERA5 already prepared and can be imported from 'ERA5/europe-era5-{year}.nc'. Please, place input cutout in the required path if not already done.")