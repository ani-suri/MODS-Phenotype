# Imports
from pathlib import Path
import jinja2
import yaml
from functools import reduce
import pyarrow as pa

run_id = "2023_07_29"

collab_path = Path("/opt/scratchspace/KLAB_SAIL/")
project_path = collab_path / "MODSPhenotypes"
mods_path = project_path / "mods"
source_path = mods_path / "src"

# YAML

## Project YAML
with open(source_path / "project_config.yaml", "rt") as fi:
    conf = fi.read().rstrip()
template = jinja2.Template(conf)
project_config = yaml.safe_load(
    template.render(
        project_path=str(project_path), run_id=str(run_id), collab_path=str(collab_path)
    )
)

## Schema YAML
with open(source_path / "encounter_schema_pandas.yaml", "rt") as fi:
    pandas_schema = yaml.safe_load(fi.read().rstrip())

with open(source_path / "encounter_schema_pyarrow.yaml", "rt") as fi:
    arrow_schema = yaml.safe_load(fi.read().rstrip())