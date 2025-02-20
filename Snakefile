
configfile: "config.yaml"

rule all:
    input: expand("sarah-era5/europe-{year}-sarah3-era5-compressed.nc", year=config["cutout"]["time"])

rule build_sara_era_cutouts:
    input: expand("sarah-era5/europe-{year}-sarah3-era5-compressed.nc", year=config["cutout"]["time"])

rule build_era5_cutouts:
    input: expand("era5/europe-era5-{year}.nc", year=config["cutout"]["time"])

rule build_sarah_cutouts:
    input: expand("sarah/europe-{year}-sarah-3.nc", year=config["cutout"]["time"])

rule build_era5_cutout:
    output: "era5/europe-era5-{year}.nc"
    threads: 2
    resources: mem_mb=32000, walltime="72:00:00" # , partition="scioi_node"
    script: "scripts/build_era5_cutout.py"

rule build_sarah_cutout:
    output: "sarah/europe-{year}-sarah-3.nc"
    threads: 2
    resources: mem_mb=32000, walltime="72:00:00" # , partition="scioi_node"
    script: "scripts/build_sarah_cutout.py"

rule merge_and_compress_cutouts:
    input:
        era5="era5/europe-era5-{year}.nc",
        sarah="sarah/europe-{year}-sarah-3.nc",
    output:
        merged="sarah-era5/europe-{year}-sarah3-era5.nc",
        merged_compressed="sarah-era5/europe-{year}-sarah3-era5-compressed.nc",
    threads: 2
    resources: mem_mb=32000, walltime="72:00:00"
    script: "scripts/merge_and_compress_cutouts.py"