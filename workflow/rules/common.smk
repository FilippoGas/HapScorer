# workflow/rules/common.smk

# import basic packages
import pandas as pd
from snakemake.utils import validate

# read sample sheet
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)["sample"]
    .tolist()
)


# validate sample sheet and config file
validate(config, schema="../schemas/config.schema.yaml")
