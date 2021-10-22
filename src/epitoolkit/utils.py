import typing as t

from pandas.core.frame import DataFrame

from .exceptions import InternalException


def extract_probes_from_region(
    manifest: DataFrame, probes: set, chr: str, start: int, end: int
) -> list:
    manifest = manifest[["CHR", "MAPINFO"]].dropna()

    probes_in_range = manifest[
        (manifest["CHR"].astype(str) == str(chr))
        & (manifest["MAPINFO"].astype(int) >= start)
        & (manifest["MAPINFO"].astype(int) <= end)
    ]
    probes_in_range = probes_in_range.sort_values("MAPINFO", ascending=True)

    if len(probes_in_range) == 0:
        raise InternalException("Probes in specific region are not available.")

    return set.intersection(set(probes_in_range.index), probes)


def check_if_porbes_in_mynorm(probes: set, cpgs_to_check: t.Collection) -> set:
    probes, cpgs_to_check = set(probes), set(cpgs_to_check)

    ovrl = set.intersection(probes, cpgs_to_check)

    if len(ovrl) != len(cpgs_to_check):
        print("Not all CpGs found in mynorm!")
        return ovrl
    else:
        return cpgs_to_check


def extract_category_from_manifest(
    manifest: DataFrame, category: str, mynorm_probes: t.Collection
) -> DataFrame:
    manifest = manifest.loc[mynorm_probes, category]
    manifest = manifest.str.split(";").explode().reset_index()

    cols = manifest.columns.tolist()
    cols[0] = "index"
    manifest.columns = cols

    manifest = manifest.drop_duplicates()
    return manifest.set_index("index")
