import typing as t
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pandas.core.frame import DataFrame

from .exceptions import InternalException


class Visualize:
    def __init__(
            self,
            manifest: t.Union[str, Path, DataFrame],
            mynorm: t.Union[str, Path, DataFrame],
            poi: t.Union[str, Path, DataFrame],
            **kwargs,
    ) -> None:

        if (
                isinstance(manifest, DataFrame)
                and isinstance(mynorm, DataFrame)
                and isinstance(poi, DataFrame)
        ):
            self.manifest = manifest
            self.mynorm = mynorm
            self.poi = poi

        elif (
                isinstance(manifest, str)
                and isinstance(mynorm, str)
                and isinstance(poi, str)
        ):
            self.manifest = pd.read_csv(
                manifest, low_memory=False, index_col=0, **kwargs
            )
            self.mynorm = pd.read_csv(mynorm, index_col=0, **kwargs)
            self.poi = pd.read_csv(poi, index_col=0, **kwargs)

        else:
            raise InternalException(
                "Both mynorm and manifest must be not empty and must be the same\
                                     type: str or DataFrame."
            )

        self.poi.columns = [column.upper() for column in self.poi.columns]
        self.probes = set.intersection(
            set(self.mynorm.index), set(self.manifest.index))

        if self.probes == {}:
            raise InternalException(
                "Any CpG overlap between myNorm and manifest files."
            )

        if not set(self.poi.index) == set(self.mynorm.columns):
            intersection = set.intersection(
                set(self.poi.index), set(self.mynorm.columns)
            )

            if len(intersection) > 0:
                print(
                    f"Not all samples in POI overlap with samples in myNorm. Future steps will use \
                       {len(intersection)} samples."
                )
                self.poi = self.poi.loc[intersection, :]
                self.mynorm = self.mynorm[intersection]

            else:
                raise InternalException(
                    "Any samples in myNorm overlap with samples in POI file."
                )

        if "POI" not in self.poi.columns:
            raise InternalException(
                "Not found columns with phenotype [POI] in POI file."
            )

    def plot_CpG(
            self,
            cpg: str,
            category_order: t.Union[t.List[str], None] = None,
            title: str = "",
            show_legend: bool = False,
            x_axis_label: str = "CpG",
            y_axis_label: str = "b-values",
            legend_title: str = "Legend",
            font_size: int = 12,
            width: int = 700,
            height: int = 500,
            y_range: t.Union[list, tuple] = (0, 1),
            static: bool = False,
            export: t.Union[str, Path, None] = None,
    ) -> None:

        if cpg not in self.probes:
            raise InternalException(f"Probes {cpg} not available.")

        data = self.mynorm.loc[cpg]
        data = pd.concat((data, self.poi), axis=1)

        fig = go.Figure()

        for sample_type in data["POI"].unique():
            values = data[data["POI"] == sample_type][cpg].values
            x_loc = [sample_type] * len(values)

            fig.add_trace(
                go.Box(
                    y=values, x=x_loc, boxpoints="all", boxmean=True, name=sample_type
                )
            )

        fig.update_layout(
            width=width,
            height=height,
            title=title,
            xaxis_title=x_axis_label,
            yaxis_title=y_axis_label,
            legend_title=legend_title,
            font=dict(size=font_size),
            showlegend=show_legend
        )

        fig.update_yaxes(range=y_range)

        if category_order:
            fig.update_xaxes(categoryorder='array', categoryarray=category_order)

        if static:
            fig.show(renderer="png")
        else:
            fig.show()

        if export and static:
            if not export.endswith(".png"):
                raise InternalException(
                    "When static export, file name must end with .png"
                )

            fig.write_image(export)

        if export and not static:
            if not export.endswith(".html"):
                raise InternalException(
                    "When interactive export, file name must end with .html"
                )

            fig.write_html(export)

    def plot_Region(self, chr: int, start: int, end: int, title: str = "",
                    show_legend: bool = False,
                    x_axis_label: str = "CpG",
                    y_axis_label: str = "b-values",
                    legend_title: str = "Legend",
                    font_size: int = 12,
                    width: int = 700,
                    height: int = 500,
                    y_range: t.Union[list, tuple] = (0, 1),
                    static: bool = False,
                    export: t.Union[str, Path, None] = None,
                    ) -> None:

        manifest = self.manifest[["CHR", "MAPINFO"]].dropna()

        probes_in_range = manifest[
            (manifest["CHR"].astype(str) == str(chr)) & (manifest["MAPINFO"].astype(int) >= start) & (
                    manifest["MAPINFO"].astype(int) <= end)]
        available_probes = set.intersection(set(probes_in_range.index), self.probes)

        if len(available_probes) == 0:
            raise InternalException("Probes in specific DMR are not available.")

        probes_in_range = probes_in_range.loc[available_probes, :].sort_values("MAPINFO", ascending=True)
        probes_in_range = list(probes_in_range.index)

        data = self.mynorm.loc[probes_in_range, :].T
        data = pd.concat((data, self.poi["POI"]), axis=1)
        data = pd.melt(data, "POI")

        fig = px.box(data, x="variable", y="value", color="POI", width=width, height=height)
        fig.update_layout(font=dict(size=font_size), xaxis_title=x_axis_label, yaxis_title=y_axis_label,
                          showlegend=show_legend, title=title, legend_title=legend_title)
        fig.update_yaxes(range=y_range)

        if static:
            fig.show(renderer="png")
        else:
            fig.show()

        if export and static:
            if not export.endswith(".png"):
                raise InternalException(
                    "When static export, file name must end with .png"
                )

            fig.write_image(export)

        if export and not static:
            if not export.endswith(".html"):
                raise InternalException(
                    "When interactive export, file name must end with .html"
                )

            fig.write_html(export)
