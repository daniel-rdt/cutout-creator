
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
        path=snakemake.output,
        module="sarah",
        sarah_dir=f"../sarah/data/sarah-3-{year}",
        chunks={"lat": -1, "time": 100},
        **cutout_params,
    )

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

    # Sarah
    if not config["existing_cutouts"]["sarah"]:
        logger.info(f"Preparing SARAH cutout for {year}")
        prepare_cutout_sarah(year, cutout_params)
        logger.info(f"SARAH cutout saved to '/sarah/europe-{year}-sarah-3.nc'")
    else:
        logger.info(f"Cutout for SARAH already prepared and can be imported from '/sarah/europe-{year}-sarah-3.nc'. Please, place input cutout in the required path if not already done.")
