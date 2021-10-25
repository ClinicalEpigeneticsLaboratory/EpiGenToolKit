import pandas as pd
from pandas.core.frame import DataFrame

import typing as t
from .exceptions import InternalException


def corect_col_names(poi: DataFrame) -> DataFrame:
    poi.columns = [column.upper() for column in poi.columns]

    if "POI" not in poi.columns:
        raise InternalException(
            "Not found columns with phenotype [POI] in POI file.")

    return poi


def find_overlap(obj_a: t.Any, obj_b: t.Any) -> set:
    try:
        return set.intersection(set(obj_a), set(obj_b))

    except:
        raise InternalException(
            "Error when trying to find overlap between objects.")


def check_objects_number(objects_collection: t.Iterable) -> None:
    if set(objects_collection) == {}:
        raise InternalException(
            "Any CpG overlap between myNorm and manifest files.")


def overlap_samples_between_poi_and_mynorm(
    poi: DataFrame, mynorm: DataFrame
) -> t.Tuple[DataFrame, DataFrame]:
    if not set(poi.index) == set(mynorm.columns):
        intersection = set.intersection(set(poi.index), set(mynorm.columns))

        if len(intersection) > 0:
            print(
                f"Not all samples in POI overlap with samples in myNorm. Future steps will use {len(intersection)} samples."
            )
            poi = poi.loc[intersection, :]
            mynorm = mynorm[intersection]
            return poi, mynorm

        else:
            raise InternalException(
                "Any samples in myNorm overlap with samples in POI file."
            )
    else:
        return poi, mynorm


def check_categories(manifest: DataFrame, user_input: t.Collection) -> None:
    if not set(user_input).issubset(set(manifest.columns)):
        diff = set.difference(set(user_input), set(manifest.columns))
        raise InternalException(f"Not found: {diff} in manifest.")
