import typing as t
from pathlib import Path

import pandas as pd
import scipy.stats as sts
import plotly.graph_objects as go
from pandas.core.frame import DataFrame

from .exceptions import InternalException
from .utils import (
    extract_probes_from_region,
    check_if_porbes_in_mynorm,
    extract_category_from_manifest,
)
from .validators import (
    corect_col_names,
    find_overlap,
    check_objects_number,
    overlap_samples_between_poi_and_mynorm,
    check_categories,
)


class Visualize:
    def __init__(
        self,
        manifest: t.Union[str, Path, DataFrame],
        mynorm: t.Union[str, Path, DataFrame],
        poi: t.Union[str, Path, DataFrame],
        poi_col: str = "POI",
        skiprows: int = 0,
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
                manifest, low_memory=False, index_col=0, skiprows=skiprows, **kwargs
            )
            self.mynorm = pd.read_csv(mynorm, index_col=0, **kwargs)
            self.poi = pd.read_csv(poi, index_col=0, **kwargs)
        else:
            raise InternalException(
                "Both mynorm and manifest must be not empty and must be the same type: str or DataFrame."
            )

        self.poi = self.poi[poi_col].to_frame()
        self.poi.columns = ["POI"]

        self.poi = corect_col_names(self.poi)
        self.probes = find_overlap(self.mynorm.index, self.manifest.index)
        check_objects_number(self.probes)

        self.poi, self.mynorm = overlap_samples_between_poi_and_mynorm(
            self.poi, self.mynorm
        )

    def __print(self, fig: go.Figure, static: bool) -> None:
        if static:
            fig.show(renderer="png")
        else:
            fig.show()

    def __upgrade_figure(
        self,
        fig: go.Figure,
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
    ) -> go.Figure:

        fig.update_layout(
            width=width,
            height=height,
            title=title,
            xaxis_title=x_axis_label,
            yaxis_title=y_axis_label,
            legend_title=legend_title,
            font=dict(size=font_size),
            showlegend=show_legend,
        )

        fig.update_yaxes(range=y_range)

        if category_order:
            fig.update_xaxes(categoryorder="array",
                             categoryarray=category_order)

        return fig

    def __export_figure(
        self,
        fig: go.Figure,
        static: bool = False,
        export: t.Union[str, Path, None] = None,
    ) -> None:
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

    def plot_CpG(
        self,
        cpg: str,
        title: str = "",
        category_order: t.Union[t.List[str], None] = None,
        static: bool = False,
        show_legend: bool = False,
        x_axis_label: str = "POI",
        y_axis_label: str = "b-values",
        legend_title: str = "Legend",
        font_size: int = 12,
        width: int = 700,
        height: int = 500,
        y_range: t.Union[list, tuple] = (0, 1),
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

        fig = self.__upgrade_figure(
            fig,
            title=title,
            show_legend=show_legend,
            x_axis_label=x_axis_label,
            y_axis_label=y_axis_label,
            legend_title=legend_title,
            font_size=font_size,
            width=width,
            height=height,
            y_range=y_range,
            category_order=category_order,
        )

        self.__print(fig, static)

        if export:
            self.__export_figure(fig, static, export)

    def plot_Range(
        self,
        chr: t.Union[str, int, None] = None,
        start: t.Union[int, None] = None,
        end: t.Union[int, None] = None,
        cpgs: t.Collection = [],
        category_order: t.Union[t.List[str], None] = None,
        show_all_points: str = "outliers",
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

        if chr and start and end:
            available_probes = extract_probes_from_region(
                self.manifest, self.probes, chr, start, end
            )

        elif cpgs:
            available_probes = check_if_porbes_in_mynorm(self.probes, cpgs)

        else:
            raise InternalException(
                "User must provide chr AND start and END arguments or any type of collection of cpgs to visualize."
            )

        available_probes = (
            self.manifest.loc[available_probes, "MAPINFO"]
            .sort_values(ascending=True)
            .index
        )
        data = self.mynorm.loc[available_probes, :].T
        data = pd.concat((data, self.poi["POI"]), axis=1)

        if category_order:
            order = category_order

        else:
            order = self.poi["POI"].unique()

        fig = go.Figure()

        for cpg in available_probes:
            for group in order:
                values = data[data["POI"] == group][cpg]
                x_loc = [cpg] * len(values)
                fig.add_trace(
                    go.Box(
                        y=values,
                        x=x_loc,
                        boxpoints=show_all_points,
                        boxmean=True,
                        name=group,
                    )
                )

        fig.update_layout(boxmode="group")

        fig = self.__upgrade_figure(
            fig,
            title=title,
            show_legend=show_legend,
            x_axis_label=x_axis_label,
            y_axis_label=y_axis_label,
            legend_title=legend_title,
            font_size=font_size,
            width=width,
            height=height,
            y_range=y_range,
        )

        self.__print(fig, static)

        if export:
            self.__export_figure(fig, static=static, export=export)


class EnrichmentAnalysis:
    def __init__(
        self,
        manifest: t.Union[str, DataFrame],
        mynorm: t.Union[str, DataFrame],
        skiprows: int = 0,
        **kwargs,
    ):

        if isinstance(manifest, DataFrame) and isinstance(mynorm, DataFrame):
            self.manifest = manifest
            self.mynorm = mynorm

        elif isinstance(manifest, str) and isinstance(mynorm, str):
            self.manifest = pd.read_csv(
                manifest, low_memory=False, index_col=0, skiprows=skiprows, **kwargs
            )
            self.mynorm = pd.read_csv(mynorm, index_col=0, **kwargs)

    @classmethod
    def load(cls, viz_object: Visualize):
        return cls(viz_object.manifest, viz_object.mynorm)

    def enrichmentAnalysis(
        self, categories_to_analyse: t.Collection, cpgs: t.Collection
    ) -> None:
        check_categories(self.manifest, categories_to_analyse)
        cpgs = check_if_porbes_in_mynorm(self.mynorm.index, cpgs)

        for category in categories_to_analyse:

            data = extract_category_from_manifest(
                self.manifest, category, self.mynorm.index
            )

            bg = data.value_counts(normalize=True, dropna=False) * 100
            input = data.loc[cpgs].value_counts(
                normalize=True, dropna=False) * 100

            frame = pd.concat((bg, input), axis=1)
            frame.columns = ["BG", "INPUT"]
            frame = frame.reset_index()
            frame = frame.set_index(category)

            _, pvalue, _, _ = sts.chi2_contingency(frame)

            print(frame)
            print(f"P-value: {pvalue}")

            fig = go.Figure(
                data=[
                    go.Bar(name="BG", x=frame.index, y=frame.BG),
                    go.Bar(name="Input", x=frame.index, y=frame.INPUT),
                ]
            )

            fig.update_layout(
                barmode="group",
                title=category,
                xaxis_title="",
                yaxis_title="Frequency [%]",
            )
            fig.show()
